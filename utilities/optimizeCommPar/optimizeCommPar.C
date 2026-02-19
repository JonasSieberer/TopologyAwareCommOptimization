/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Generate ranklist to optimize rank to core assignment according
    to communication costs at the processor patches and the given
    NUMA domains in the decompositionParDict.
    Create optimized communication schedule for global reduce methods.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvMesh.H"
#include "processorRunTimes.H"
#include "topoSetSource.H"
#include "volFields.H"
#include "systemDict.H"
#include "decompositionMethod.H"
#include "processorFvPatch.H"
#include "CompactListList.H"
#include "polyBoundaryMeshEntries.H"
#include "scotch.H"
#include "OFstream.H"
#include "SortableList.H"

#include <asm-generic/errno.h>
#include <cwctype>
#include <fstream>
#include <cstdio>
#include <pthread.h>
#include <regex>
#include <algorithm>
#include <cstdio>   // For popen()
#include <cstdlib>  // For atoi()

#include "IFstream.H"
#include "IStringStream.H"

#include <iostream>
#include <chrono>
#include <unistd.h>
#include <sstream>
#include <iomanip>

extern "C"
{
#include "scotch.h"
}

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Configuration structure to use in whole utility
struct NumaConfig
{
    label nDomains;         // number of subdomains
    label nodes;            // number of used nodes
    label procPN;           // number of cores per node
    label numaNodeSizes;    // size of a NUMA node
    label numaDomains;      // number of NUMA domains
    label numaNodes;        // number of NUMA nodes

    label version;          // only global or local+global version

    dictionary decomposeParDict;
    dictionary optimizeCoeffsDict;

    labelListList numa;     //
    labelList numaList;     // read in NUMA nodes
    labelList domainList;   // read in NUMA domains
    labelList result;       // rank to node assignment
    labelListList rankList; // final optimized ranList to print as rankFile

    // Default constructor
    NumaConfig()
    :
        nDomains(0),
        nodes(0),
        procPN(0),
        numaNodeSizes(0),
        numaDomains(0),
        numaNodes(0),
        version(0),
        decomposeParDict(),
        optimizeCoeffsDict(),
        numa(),
        numaList(),
        domainList(),
        result(),
        rankList()
    {}
};

NumaConfig cfg;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Error-handling for decompose method with Scotch
void checkScotch
(
    const int retVal,
    const char* str
)
{
    if (retVal)
    {
        FatalErrorInFunction
            << "Call to scotch routine " << str << " failed."
            << exit(FatalError);
    }
}

// Helper method for balancing partitons:
label exchange
(
    const labelList& adjncy,
    const labelList& xadj,
    const scalarField& edgeWeights,

    labelList& subSizes,
    labelList& decomp,

    const label& from,
    const label& to,
    const label& weight
)
{
    label minWeight = 0;
    bool first = true;
    label proc = -1;

    // Iterate over all possible processes on "from" (edge to "to")
    // and check which leads to minimal new edge weight
    for (label i = 0; i < decomp.size(); i++)
    {
        if(decomp[i] == from)
        {


            label nbrWeight = 0;
            label ownWeight = 0;

            for (label j = xadj[i]; j < xadj[i+1]; j++)
            {
                if (decomp[adjncy[j]] == to)
                {
                    nbrWeight += edgeWeights[j];
                }

                if (decomp[adjncy[j]] == from)
                {
                    ownWeight += edgeWeights[j];
                }
            }

            if (nbrWeight > 0)
            {
                if (first)
                {
                    minWeight = ownWeight - nbrWeight;
                    proc = i;
                    first = false;
                }
                else
                {
                    if (ownWeight - nbrWeight < minWeight)
                    {
                        minWeight = ownWeight - nbrWeight;
                        proc = i;
                    }
                }
            }
        }
    }

    // Exchange proc and update partition

    Info << "\tMove proc " << proc << " from " << from << " to " << to << endl;

    label newWeight = weight + minWeight;

    decomp[proc] = to;
    subSizes[from]--;
    subSizes[to]++;

    return newWeight;
}

// Simple graph balancer for not balanced partitions
label balanceDecomp
(
    const labelList& adjncy,
    const labelList& xadj,
    const scalarField& edgeWeights,

    const label& maxSize,
    const label& parts,
    const labelList& partSizes,
    labelList& subSizes,

    labelList& decomp
)
{
    Info << "\nBalance Partition:\n\n";

    Info << "Current partition sizes: " << subSizes << endl;

    // Create graph of the partition with each part a node and edges showing adjancies
    // Calculate all surplus and defencies to calculate flows

    labelList surpDef = subSizes - partSizes;

    labelListList nbrPartWeights(parts, labelList(parts, 0));

    for (label verti = 0; verti < decomp.size(); verti++)
    {
        for (label i = xadj[verti]; i < xadj[verti+1]; i++)
        {
            nbrPartWeights[decomp[verti]][decomp[adjncy[i]]] += edgeWeights[i];
        }
    }

    labelList qXadj(parts+1, 0);

    for (label i = 0; i < parts; i++)
    {
        label nbrs = 0;
        for (label j = 0; j < parts; j++)
        {
            if (i != j)
            {
                if (nbrPartWeights[i][j] > 0)
                {
                    nbrs++;
                }
            }
        }

        qXadj[i+1] = qXadj[i] + nbrs;
    }

    labelList qAdjncy(qXadj[parts]);
    labelList qEdgeWeights(qXadj[parts]);

    label count = 0;
    for (label i = 0; i < parts; i++)
    {
        for (label j = 0; j < parts; j++)
        {
            if (i != j)
            {
                if (nbrPartWeights[i][j] > 0)
                {
                    qAdjncy[count] = j;
                    qEdgeWeights[count] = nbrPartWeights[i][j];
                    count++;
                }
            }
        }
    }

    bool balanced = false;

    // While it is not balanced account one surplus after the other

    while (!balanced)
    {
        label surp = -1;

        for (label i = 0; i < parts; i++)
        {
            if (surpDef[i] > 0)
            {
                surp = i;
            }
        }

        // Perform a breadth first search on the whole graph

        labelList parent(parts, -1);
        labelList queue(parts, -1);
        label head = 0;
        label tail = 0;

        parent[surp] = surp;
        queue[tail] = surp;
        tail++;

        while (head < tail)
        {
            label curr = queue[head];

            for (label i = qXadj[curr]; i < qXadj[curr+1]; i++)
            {
                if (parent[qAdjncy[i]] == -1)
                {
                    parent[qAdjncy[i]] = curr;
                    queue[tail] = qAdjncy[i];
                    tail++;
                }
            }

            head++;
        }

        labelList length(parts, 0);

        for (label i = 0; i < parts; i++)
        {
            if (i == surp)
            {
                continue;
            }

            label path = parent[i];
            length[i]++;

            while (path != surp)
            {
                path = parent[path];
                length[i]++;
            }
        }

        // Find nearest deficit

        label minLength = parts;
        label def = -1;

        for (label i = 0; i < parts; i++)
        {
            if (surpDef[i] < 0)
            {
                if (length[i] < minLength)
                {
                    def = i;
                    minLength = length[i];
                }
            }
        }

        // Create path from surplus to deficit

        labelList path(length[def]+1, -1);
        path[0] = def;

        for (label i = 1; i < path.size(); i++)
        {
            path[i] = parent[path[i-1]];
        }

        // Exchange procs along the path with exchange method
        // which determines the optimal proc to exchange

        Info << "Exchange one vertex on the path: ";
        for (label i = length[def]; i > 0; i--)
        {
            Info << path[i] << " -> ";
        }
        Info << path[0] << endl;


        for (label i = length[def]; i > 0; i--)
        {
            label from = path[i];
            label to = path[i-1];

            label newWeight = exchange(
                adjncy,
                xadj,
                edgeWeights,
                subSizes,
                decomp,
                from,
                to,
                nbrPartWeights[from][to]
            );

            nbrPartWeights[from][to] = newWeight;
            nbrPartWeights[to][from] = newWeight;

            surpDef = subSizes - partSizes;

        }

        // Adapt partition size and check if it is balanced now

        balanced = true;

        for (label i = 0; i < parts; i++)
        {
            if (subSizes[i] != maxSize)
            {
                balanced = false;
            }
        }
    }

    Info << "Partition now balanced!\n\n";

    return 0;
}


// decompose graph with Scotch
// copied and adapted from OpenFOAM
label decomposeScotch
(
    const dictionary& dict,
    const fileName& meshPath,
    const labelList& adjncy,
    const labelList& xadj,
    const scalarField& edgeWeights,

    const label& maxSize,
    const label& parts,
    const labelList& partSizes,

    labelList& decomp
)
{
    // Dump graph
    if (dict.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            dict.subDict("scotchCoeffs");

        if (scotchCoeffs.lookupOrDefault("writeGraph", false))
        {
            OFstream str(meshPath + ".grf");

            Info<< "Dumping Scotch graph file to " << str.name() << endl
                << "Use this in combination with gpart." << endl;

            label version = 0;
            str << version << nl;
            // Number of vertices
            str << xadj.size()-1 << ' ' << adjncy.size() << nl;
            // Numbering starts from 0
            label baseval = 0;
            // Has weights?
            label hasEdgeWeights = 1;
            label hasVertexWeights = 0;
            label numericflag = 10*hasEdgeWeights+hasVertexWeights;
            str << baseval << ' ' << numericflag << nl;
            for (label celli = 0; celli < xadj.size()-1; celli++)
            {
                label start = xadj[celli];
                label end = xadj[celli+1];
                str << end-start;

                for (label i = start; i < end; i++)
                {
                    str << ' ' << adjncy[i];
                }
                str << nl;
            }
        }
    }

    // Reset the seed of the pseudo-random generator used by the graph
    // partitioning routines of the libScotch library. Two consecutive calls to
    // the same libScotch partitioning routines, and separated by a call to
    // SCOTCH randomReset, will always yield the same results, as if the
    // equivalent standalone Scotch programs were used twice, independently,
    SCOTCH_randomReset();

    // Strategy
    // ~~~~~~~~

    // Default.
    SCOTCH_Strat stradat;
    checkScotch(SCOTCH_stratInit(&stradat), "SCOTCH_stratInit");

    if (dict.found("scotchCoeffs"))
    {
        const dictionary& scotchCoeffs =
            dict.subDict("scotchCoeffs");

        string strategy;
        if (scotchCoeffs.readIfPresent("strategy", strategy))
        {
            SCOTCH_stratGraphMap(&stradat, strategy.c_str());
        }
    }

    // Graph
    // ~~~~~

    // Graph Initialization
    SCOTCH_Graph graph;
    checkScotch(SCOTCH_graphInit(&graph), "SCOTCH_graphInit");

    labelList edlotab;
    for (label i = 0; i < edgeWeights.size(); i++)
    {
        edlotab.append(edgeWeights[i]);
    }

    labelList vendtab;
    for (label i = 1; i < xadj.size(); i++)
    {
        vendtab.append(xadj[i]);
    }

    SCOTCH_Graph grafdat;
    checkScotch(SCOTCH_graphInit(&grafdat), "SCOTCH_graphInit");
    checkScotch
    (
        SCOTCH_graphBuild
        (
            &grafdat,
            0,                      // baseval, c-style numbering
            xadj.size()-1,          // vertnbr, nCells
            xadj.begin(),           // verttab, start index per cell into adjncy
            vendtab.begin(),        // vendtab, end index  ,,
            nullptr,                // velotab, vertex weights
            nullptr,                // vlbltab
            adjncy.size(),          // edgenbr, number of arcs
            adjncy.begin(),         // edgetab
            edlotab.begin()         // edlotab, edge weights
        ),
        "SCOTCH_graphBuild"
    );
    checkScotch(SCOTCH_graphCheck(&grafdat), "SCOTCH_graphCheck");

    // Architecture
    SCOTCH_Arch archdat;
    checkScotch(SCOTCH_archInit(&archdat), "SCOTCH_archInit");

    decomp.setSize(xadj.size()-1);
    decomp = 0;

    checkScotch
    (
        SCOTCH_archCmpltw(&archdat, parts, partSizes.begin()),
        "SCOTCH_archCmpltw"
    );

    checkScotch(
        SCOTCH_graphMap(
            &grafdat,
            &archdat,
            &stradat,
            decomp.begin()
        ),
        "SCOTCH_graphMap"
    );

    labelList subSizes(parts, 0);
    for (label i = 0; i < decomp.size(); i++)
    {
        subSizes[decomp[i]]++;
    }

    // Balance decomp if part sizes are not equal
    for (label parti = 0; parti < parts; parti++)
    {
        if (subSizes[parti] > maxSize)
        {
           balanceDecomp(
                adjncy,
                xadj,
                edgeWeights,
                maxSize,
                parts,
                partSizes,
                subSizes,
                decomp
            );
        }
    }

    labelList sizes(parts, 0);
    for (label i = 0; i < decomp.size(); i++)
    {
        sizes[decomp[i]]++;
    }

    #ifdef FE_NOMASK_ENV
    feenableexcept(oldExcepts);
    #endif

    // Release storage for graph
    SCOTCH_graphExit(&grafdat);

    // Release storage for strategy
    SCOTCH_stratExit(&stradat);

    // Release storage for network topology
    SCOTCH_archExit(&archdat);

    return 0;
}

// Read in the node names from the dictionary
// Or get node names from SLURM command
// For the rankFile the exact names of the used nodes are important
List<string> readNodeList(const dictionary& dict)
{
    List<string> nodeList = dict.lookupOrDefault<List<string>>("nodeList", List<string>());

    if (!nodeList.empty())
    {
        return nodeList;
    }

    const char* slurmNodeList = std::getenv("SLURM_NODELIST");

    if (!slurmNodeList) {
        Foam::FatalError << "No nodeList could be found.\n"
                << "Provide nodeList in allocateParDict or"
                << " via a job script. "
                << Foam::abort(Foam::FatalError);
    }

    // Create the scontrol command to expand node names
    string command = "scontrol show hostnames " + std::string(slurmNodeList);
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe) {
        Foam::FatalError << "Failed to run scontrol command.";
    }

    char buffer[128];
    while (fgets(buffer, sizeof(buffer), pipe.get()) != nullptr) {
        string node(buffer);
        node.erase(node.find_last_not_of(" \n\r\t") + 1); // Remove trailing newline
        nodeList.append(node);
    }

    return nodeList;
}


// Read in decomposition graph from patch description in ./constant/polyMesh/boundary
void readDecomp(
    labelList& adjncy,
    labelList& xadj,
    scalarField& edgeWeights,
    const Time& runTime
)
{
    const fileName caseDir = getEnv("FOAM_CASE");
    const fileName caseName = getEnv("FOAM_CASENAME");

    label edgeCount = 0;

    xadj.append(0);

    for (label proci = 0; proci < cfg.nDomains; proci++)
    {
        fileName meshDir(caseDir + "/processor" + name(proci));
        fileName boundaryPath = meshDir/"constant/polyMesh";

        if (!isFile(boundaryPath/"boundary"))
        {
            FatalErrorInFunction
                << "Boundary file does not exist: " << boundaryPath << nl
                << "Check if decomposePar was run correctly." << exit(FatalError);
        }

        typeIOobject<polyBoundaryMesh> ioObj
        (
            "boundary",
            boundaryPath,
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );


        if (ioObj.headerOk())
        {
            polyBoundaryMeshEntries patchEntries(ioObj);

            forAll(patchEntries, patchi)
            {
                word type(patchEntries[patchi].dict().lookup<word>("type"));

                if (type == "processor")
                {
                    adjncy.append(patchEntries[patchi].dict().lookup<label>("neighbProcNo"));
                    edgeWeights.append(patchEntries[patchi].dict().lookup<label>("nFaces"));

                    edgeCount++;
                }
            }

            xadj.append(edgeCount);
        }
        else
        {
            FatalErrorInFunction
                << "Cannot find or read the boundary file for case "
                << nl
                << "    " << meshDir
                << exit(FatalError);
        }
    }

    return;
}

// Calculate amount of processes for each node
void fillAmount(labelList& partSizes, label totKnots, label& maxPerPart)
{
    label count = 0;
    label restDom = totKnots;

    while (restDom > maxPerPart)
    {
        partSizes[count] = maxPerPart;
        count++;
        restDom -= maxPerPart;
    }

    partSizes[count] = restDom;
    partSizes.resize(count+1);
}

// Creates the collected rankList from the single assignment lists
void createRankList(labelList& result, labelList& domainResult, labelList& numaResult)
{
    labelListList checkNuma(cfg.nodes, labelList(cfg.procPN, -1));
    for (label i = 0; i < cfg.numa.size(); i++)
    {
        for (label j = 0; j < cfg.numa[i].size(); j++)
        {
            checkNuma[i][j] = cfg.numa[i][j];
        }
    }

    labelListList checkDomain(cfg.nodes, labelList(cfg.procPN, -1));
    for (label i = 0; i < cfg.nodes; i++)
    {
        for (label j = 0; j < cfg.domainList.size(); j++)
        {
            checkDomain[i][j] = cfg.domainList[j];
        }
    }

    for (label i = 0; i < cfg.nDomains; i++)
    {
        label node = result[i];
        label domain = domainResult[i];
        label numa = numaResult[i];

        for (label j = 0; j < cfg.procPN; j++)
        {
            if (checkDomain[node][j] == domain)
            {
                if (checkNuma[node][j] == numa)
                {
                    numa = j;
                    checkDomain[node][j] = -1;
                    checkNuma[node][j] = -1;
                    break;
                }
            }
        }

        cfg.rankList[0][i] = node;
        cfg.rankList[1][i] = numa;
    }
}

// Sort the results such that rank 0 is assigned to the first slot on the first node
// and the following hardware hierachies are assigned ascending, without altering the
// described communication structure
void sortResults(labelList& result, labelList& domainResult, labelList& numaResult)
{
    label cnode = 0;
    for (label i = 0; i < cfg.nDomains; i++)
    {
        if (result[i] > cnode)
        {
            label a = result[i];
            label b = cnode;

            for (label j = 0; j < cfg.nDomains; j++)
            {
                if (result[j] == a)
                {
                    result[j] = b;
                }
                else if (result[j] == b)
                {
                    result[j] = a;
                }
            }

            cnode++;
        }

        if (result[i] == cnode)
        {
            cnode++;
        }
    }

    labelList cdomain(cfg.nodes, 0);
    for (label i = 0; i < cfg.nDomains; i++)
    {
        if (domainResult[i] > cdomain[result[i]])
        {
            label a = domainResult[i];
            label b = cdomain[result[i]];

            for (label j = 0; j < cfg.nDomains; j++)
            {
                if (domainResult[j] == a && result[j] == result[i])
                {
                    domainResult[j] = b;
                }
                else if (domainResult[j] == b && result[j] == result[i])
                {
                    domainResult[j] = a;
                }
            }

            cdomain[result[i]]++;
        }

        if (domainResult[i] == cdomain[result[i]])
        {
            cdomain[result[i]]++;
        }
    }

    labelListList cnuma(cfg.nodes, labelList(cfg.numaDomains, 0));
    for (label i = 0; i < cfg.nDomains; i++)
    {
        if (numaResult[i] > cnuma[result[i]][domainResult[i]])
        {
            label a = numaResult[i];
            label b = cnuma[result[i]][domainResult[i]];

            for (label j = 0; j < cfg.nDomains; j++)
            {
                if (numaResult[j] == a && result[j] == result[i] && domainResult[i] == domainResult[j])
                {
                    numaResult[j] = b;
                }
                else if (numaResult[j] == b && result[j] == result[i] && domainResult[i] == domainResult[j])
                {
                    numaResult[j] = a;
                }
            }

            cnuma[result[i]][domainResult[i]]++;
        }

        if (numaResult[i] == cnuma[result[i]][domainResult[i]])
        {
            cnuma[result[i]][domainResult[i]]++;
        }
    }

    return;
}

// Printing the rankList as ./constant/rankFile
void printRankFile
(
    const List<string>& nodeList,
    const Time& runTime
)
{
    fileName rankFilePath(runTime.constant()/"rankFile");

    Info<< "Writing MPI rankfile to " << rankFilePath << nl;

    OFstream rankFile(rankFilePath);

    if (!rankFile.good())
    {
        FatalErrorInFunction
            << "Failed to create rankfile at "
            << rankFilePath << exit(FatalError);
    }

    for (label i = 0; i < cfg.nDomains; i++)
    {
        rankFile << "rank " << i << "=" << nodeList[cfg.rankList[0][i]].c_str();
        rankFile << " slot=" << cfg.rankList[1][i] << "\n";
    }
}

// Helper method to get rank and core from read in line string
bool parseRankAndCore(const Foam::string& line, label& rank, label& core)
{
    using sizeType = Foam::string::size_type;

    sizeType posRank = line.find("rank ");
    if (posRank == Foam::string::npos)
    {
        return false;
    }

    posRank += 5; // skip "rank "
    sizeType endRank = line.find(" ", posRank);
    if (endRank == Foam::string::npos)
    {
        return false;
    }

    {
        IStringStream is(line.substr(posRank, endRank - posRank));
        is >> rank;
    }

    sizeType posCore = line.find("core ");
    if (posCore == Foam::string::npos)
    {
        return false;
    }

    posCore += 5; // skip "core "
    sizeType endCore = line.find("[", posCore);
    if (endCore == Foam::string::npos)
    {
        return false;
    }

    {
        IStringStream is(line.substr(posCore, endCore - posCore));
        is >> core;
    }

    return true;
}

// For generating only the optimized global communication schedule for OpenFOAM
// the rank to core assignment from the present used MPI version has to be read in
// from a file previous generated by a test run with the flag "--report-bindings"
void readBindingFile(const fileName& file, labelList& rankToCore)
{
    IFstream is(file);
    if (!is.good())
    {
        FatalErrorInFunction << "Cannot open file: " << file << exit(FatalError);
    }

    string line;

    labelList corr(rankToCore.size(), 0);

    while (is.good())
    {
        is.getLine(line);
        if (line.empty()) continue;

        label rank = -1;
        label core = -1;

        if (parseRankAndCore(line, rank, core))
        {
            if (rank < cfg.procPN)
            {
                rankToCore[rank] = core;
                corr[rank] = 1;
            }
        }
    }

    label sum = 0;
    for (label i = 0; i < corr.size(); i++)
    {
        sum += corr[i];
    }

    if (sum != cfg.procPN)
    {
        FatalErrorInFunction
                << "Wrong reportBindings.txt file!" << nl
                << exit(FatalError);
    }

    return;
}

// Calculate the optimized global communication strucure.
// Create a total schedule by defining the OpenFOAM tree structure
// on each hierachical layer.
void calcCommunicator
(
    labelList& result,
    labelList& domainResult,
    labelList& numaResult,
    const Time& runTime
)
{
    // For statistics count the communication on all comm types for global comm
    // (inter node, inter numa domain, inter numa node, intra numa node)
    labelList commCount(4, 0);

    // Map results on common nodes, numa domains and numa nodes
    List<List<List<DynamicList<label>>>> mapping(cfg.nodes, List<List<DynamicList<label>>>
                                                 (cfg.numaDomains, List<DynamicList<label>>
                                                 (cfg.numaNodes, DynamicList<label>()
                                                )));

    for (label ranki = 0; ranki < cfg.nDomains; ranki++)
    {
        mapping[result[ranki]][domainResult[ranki]][numaResult[ranki]].append(ranki);
    }

    // Build tree communication structure on all three common communication types
    List<DynamicList<label>> receives(cfg.nDomains);
    labelList sends(cfg.nDomains, -1);


    for (label i = 0; i < mapping.size(); i++)
    {
        for (label j = 0; j < mapping[i].size(); j++)
        {
            for (label k = 0; k < mapping[i][j].size(); k++)
            {
                // tree communication on all numa nodes
                label thisNumaNodeSize = mapping[i][j][k].size();
                label nLevels = 1;
                while ((1 << nLevels) < thisNumaNodeSize)
                {
                    nLevels++;
                }

                label offset = 2;
                label childOffset = offset/2;

                for (label level = 0; level < nLevels; level++)
                {
                    label receiveID = 0;
                    while (receiveID < thisNumaNodeSize)
                    {
                        // Determine processor that sends and we receive from
                        label sendID = receiveID + childOffset;

                        if (sendID < thisNumaNodeSize)
                        {
                            receives[mapping[i][j][k][receiveID]].append(mapping[i][j][k][sendID]);
                            sends[mapping[i][j][k][sendID]] = mapping[i][j][k][receiveID];

                            commCount[3]++;
                        }

                        receiveID += offset;
                    }

                    offset <<= 1;
                    childOffset <<= 1;
                }
            }

            // tree communication on all numa domains
            label thisNumaDomainSize = mapping[i][j].size();
            label nLevels = 1;
            while ((1 << nLevels) < thisNumaDomainSize)
            {
                nLevels++;
            }

            label offset = 2;
            label childOffset = offset/2;

            for (label level = 0; level < nLevels; level++)
            {
                label receiveID = 0;
                while (receiveID < thisNumaDomainSize)
                {
                    // Determine processor that sends and we receive from
                    label sendID = receiveID + childOffset;

                    if (sendID < thisNumaDomainSize)
                    {
                        receives[mapping[i][j][receiveID][0]].append(mapping[i][j][sendID][0]);
                        sends[mapping[i][j][sendID][0]] = mapping[i][j][receiveID][0];

                        commCount[2]++;
                    }

                    receiveID += offset;
                }

                offset <<= 1;
                childOffset <<= 1;
            }
        }

        // tree communication on all nodes
        label thisNodeSize = mapping[i].size();
        label nLevels = 1;
        while ((1 << nLevels) < thisNodeSize)
        {
            nLevels++;
        }

        label offset = 2;
        label childOffset = offset/2;

        for (label level = 0; level < nLevels; level++)
        {
            label receiveID = 0;
            while (receiveID < thisNodeSize)
            {
                // Determine processor that sends and we receive from
                label sendID = receiveID + childOffset;

                if (sendID < thisNodeSize)
                {
                    receives[mapping[i][receiveID][0][0]].append(mapping[i][sendID][0][0]);
                    sends[mapping[i][sendID][0][0]] = mapping[i][receiveID][0][0];

                    commCount[1]++;
                }

                receiveID += offset;
            }

            offset <<= 1;
            childOffset <<= 1;
        }

    }

    // tree communication on root node
    label nLevels = 1;
    while ((1 << nLevels) < cfg.nodes)
    {
        nLevels++;
    }


    label offset = 2;
    label childOffset = offset/2;

    for (label level = 0; level < nLevels; level++)
    {
        label receiveID = 0;
        while (receiveID < cfg.nodes)
        {
            // Determine processor that sends and we receive from
            label sendID = receiveID + childOffset;

            if (sendID < cfg.nodes)
            {
                receives[mapping[receiveID][0][0][0]].append(mapping[sendID][0][0][0]);
                sends[mapping[sendID][0][0][0]] = mapping[receiveID][0][0][0];

                commCount[0]++;
            }

            receiveID += offset;
        }

        offset <<= 1;
        childOffset <<= 1;
    }

    // write the communication shedula as a dictionary into constant
    Info<< "Writing optimized communication schedule to \"constant/communicationSchedule\"" << nl;

    IOdictionary dict
    (
        IOobject
        (
            "communicationSchedule",
            runTime.constant(),      // <case>/constant
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    dict.add("nDomains", sends.size());
    dict.add("sends", sends);
    dict.add("receives", receives);

    dict.regIOobject::write();


    // Replicate the OpenFOAM standard tree communication schedule to get comparabel statistics
    labelList commCountStandard(4, 0);
    {
        label nProcs = cfg.nDomains;
        label nLevels = 1;
        while ((1 << nLevels) < nProcs)
        {
            nLevels++;
        }

        label offset = 2;
        label childOffset = offset/2;

        for (label level = 0; level < nLevels; level++)
        {
            label receiveID = 0;
            while (receiveID < nProcs)
            {
                // Determine processor that sends and we receive from
                label sendID = receiveID + childOffset;

                if (sendID < nProcs)
                {
                    if (cfg.result[sendID] != cfg.result[receiveID])
                    {
                        commCountStandard[0]++;
                    }
                    else
                    {
                        if (domainResult[sendID] != domainResult[receiveID])
                        {
                            commCountStandard[1]++;
                        }
                        else
                        {
                            if (numaResult[sendID] != numaResult[receiveID])
                            {
                                commCountStandard[2]++;
                            }
                            else
                            {
                                commCountStandard[3]++;
                            }
                        }
                    }
                }

                receiveID += offset;
            }

            offset <<= 1;
            childOffset <<= 1;
        }
    }

    Info << "\nglobal communication statistics:\n";
    Info << "(inter node, inter numa domain, inter numa node, intra numa node)\n";
    Info << "optimized schedule: (" << commCount[0] << ", " << commCount[1] << ", ";
    Info << commCount[2] << ", " << commCount[3] << ")\n";
    Info << "std tree schedule:  (" << commCountStandard[0] << ", " << commCountStandard[1] << ", ";
    Info << commCountStandard[2] << ", " << commCountStandard[3] << ")\n\n";

    return;
}

// Calculate statistics for local communication and print to log file
void calcStatistics(
    labelList& domainResult,
    labelList& numaResult,
    labelList& adjncy,
    labelList& xadj,
    scalarField& edgeWeights
)
{
    // Count amount of edges and their weights for all comm types of local comm
    // (inter node, inter numa domain, inter numa node, intra numa node)
    labelList edgeCount(4, 0);
    scalarList edgeWeight(4, 0.0);

    for (label send = 0; send < xadj.size()-1; send++)
    {
        for (label recvi = xadj[send]; recvi < xadj[send+1]; recvi++)
        {
            label recv = adjncy[recvi];

            if (cfg.result[send] != cfg.result[recv])
            {
                edgeCount[0]++;
                edgeWeight[0] += edgeWeights[recvi];
            }
            else
            {
                if (domainResult[send] != domainResult[recv])
                {
                    edgeCount[1]++;
                    edgeWeight[1] += edgeWeights[recvi];
                }
                else
                {
                    if (numaResult[send] != numaResult[recv])
                    {
                        edgeCount[2]++;
                        edgeWeight[2] += edgeWeights[recvi];
                    }
                    else
                    {
                        edgeCount[3]++;
                        edgeWeight[3] += edgeWeights[recvi];
                    }
                }
            }
        }
    }

    for (label i = 0; i < 4; i++)
    {
        edgeCount[i] /= 2;
        edgeWeight[i] /=2;
    }

    Info << "\nlocal communication statistics:\n";
    Info << "(inter node, inter numa domain, inter numa node, intra numa node)\n";
    Info << "edgeCount:  (" << edgeCount[0] << ", " << edgeCount[1] << ", ";
    Info << edgeCount[2] << ", " << edgeCount[3] << ")\n";
    Info << "edgeWeight: (" << edgeWeight[0] << ", " << edgeWeight[1] << ", ";
    Info << edgeWeight[2] << ", " << edgeWeight[3] << ")\n\n";

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Initialise utility and read in numa configuration
    #include "setRootCase.H"

    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "optimizeCommPar must be run with -parallel"
            << exit(FatalError);
    }

    // Set time from database
    processorRunTimes runTimes(Foam::Time::controlDictName, args);
    const Time& runTime = runTimes.completeTime();

    // Get the decomposition dictionary
    cfg.decomposeParDict = decompositionMethod::decomposeParDict(runTimes.completeTime());

    if (!cfg.decomposeParDict.found("optimizeCommCoeffs"))
    {
        FatalErrorInFunction
                << "No coefficients to execute optimizeCommPar could be found" << nl
                << exit(FatalError);
    }

    // Get number of Subdomains
    cfg.nDomains = cfg.decomposeParDict.lookup<label>("numberOfSubdomains");

    // Get the subdictionary 'optimizeCommCoeffs'
    cfg.optimizeCoeffsDict = cfg.decomposeParDict.subDict("optimizeCommCoeffs");

    // Get the node number list from the subdictionary
    List<string> nodeList(readNodeList(cfg.optimizeCoeffsDict));

    cfg.nodes = nodeList.size();

    // Get the nodes from the subdictionary
    cfg.version = cfg.optimizeCoeffsDict.lookupOrDefault<label>("version", 0);

    // Get the architecture of each node from subdirectory
    cfg.procPN = cfg.optimizeCoeffsDict.lookup<label>("procPerNode");

    // Get the number of numa domains from the subdictionary
    cfg.numaDomains = cfg.optimizeCoeffsDict.lookup<label>("numaDomains");

    // Get the number of numa nodes from the subdictionary
    cfg.numaNodes = cfg.optimizeCoeffsDict.lookup<label>("numaNodes");

    cfg.numaNodeSizes = cfg.procPN / cfg.numaDomains / cfg.numaNodes;

    if (cfg.version > 2)
    {
        FatalErrorInFunction
                << "Not a correct version number. Possible versions:" << nl
                << "0: Use scotch to create rankfile from read in numa nodes and domains" << nl
                << "1: Random rankfile allocation" << nl
                << "2: Only create communication schedule on OpenFOAM standard with read in numa nodes and domains"
                << nl
                << exit(FatalError);
    }

    Info << "nDomains: " << cfg.nDomains << endl;
    Info << "version: ";
    if (cfg.version == 0) {Info << "Scotch on read in numa nodes/domains" << endl;}
    else if (cfg.version == 1) {Info << "Random allocation" << endl;}
    else if (cfg.version == 2) {Info << "Only communication schedule on OpenFOAM standard";
        Info << " from read in numa nodes/domains" << endl;}
    else {Info << "Wrong version coefficient" << endl;}

    Info << "nodes: " << cfg.nodes << endl;
    Info << "processes per Node: " << cfg.procPN << endl;
    Info << "numaDomains: " << cfg.numaDomains << endl;
    Info << "numaNodes: " << cfg.numaNodes << endl;
    Info << "numaNodeSizes: " << cfg.numaNodeSizes << endl;
    Info << "nodeList: " << nodeList << endl;

    // Check if enough cores are used
    if (cfg.nodes * cfg.procPN < cfg.nDomains)
    {
        FatalErrorInFunction
                << "Not enough cores and/or nodes " << nl
                << "Running " << cfg.nDomains << " processes on " << cfg.nodes
                << " nodes with " << cfg.procPN << " cores per node" << nl
                << "Max possible processes here are " << cfg.nodes*cfg.procPN << nl
                << "Use more nodes or run with fewer processes"
                << nl
                << exit(FatalError);
    }

    // Check the specified number of processes is consistent with any existing
    // processor directories
    {
        const label nProcs0 =
            fileHandler().nProcs(runTimes.completeTime().path());

        if (nProcs0 && nProcs0 != cfg.nDomains)
        {
            FatalErrorInFunction
                << "Case is already decomposed with " << nProcs0
                << " domains, instead of " << cfg.nDomains
                << " domains as specified in decomposeParDict"
                << nl
                << exit(FatalError);
        }
    }

    // Read in the decomposition as a graph
    labelList adjncy, xadj;
    scalarField edgeWeights;

    // Read in decomposition and save as graph
    readDecomp(adjncy, xadj, edgeWeights, runTime);

    cfg.result = labelList(cfg.nDomains, 0);
    cfg.rankList = labelListList(2, labelList(cfg.nDomains, -1));

    if (Pstream::myProcNo() == 0)
    {

    labelList domainResult(cfg.nDomains, -1);
    labelList numaResult(cfg.nDomains, -1);

    // Read in numaList and domainList (for all versions -> calc statistics)
    {
        // Get the numa list from the subdictionary
        cfg.numaList = cfg.optimizeCoeffsDict.lookupOrDefault<labelList>("numaList", labelList());

        if (cfg.numaList.size() == 0)
        {
            FatalErrorInFunction
                    << "No descripted numa list found in decomposeParDict " << nl
                    << exit(FatalError);
        }

        if (cfg.numaList.size() != cfg.procPN)
        {
            FatalErrorInFunction
                    << "Wrong descripted numa list found in decomposeParDict " << nl
                    << "List describes " << cfg.numaList.size() << " processes per node, instead of"
                    << cfg.procPN << " processes per node" << nl
                    << exit(FatalError);
        }

        cfg.numa = labelListList(cfg.nodes, labelList(cfg.procPN, 0));

        for (label i = 0; i < cfg.numa.size(); i++)
        {
            for (label j = 0; j < cfg.procPN; j++)
            {
                cfg.numa[i][j] = cfg.numaList[j];
            }
        }

        // Get the numa list from the subdictionary
        cfg.domainList = cfg.optimizeCoeffsDict.lookupOrDefault<labelList>("domainList", labelList());

        if (cfg.domainList.size() == 0 && cfg.numaDomains > 1)
        {
            FatalErrorInFunction
                    << "No descripted domain list found in decomposeParDict " << nl
                    << exit(FatalError);
        }

        if (cfg.domainList.size() != cfg.procPN)
        {
            FatalErrorInFunction
                    << "Wrong descripted domain list found in decomposeParDict " << nl
                    << "List describes " << cfg.numaList.size() << " processes per node, instead of"
                    << cfg.procPN << " processes per node" << nl
                    << exit(FatalError);
        }
    }

    // Write rank file for random allocation
    if (cfg.version == 1)
    {
        labelList random(cfg.nDomains);
        for (label i = 0; i < cfg.nDomains; i++)
        {
            random[i] = i;
        }
        std::random_shuffle(std::begin(random), std::end(random));

        for (label i = 0; i < cfg.nDomains; i++)
        {
            label node = 0;

            while (random[i] >= cfg.procPN)
            {
                random[i] -= cfg.procPN;
                node++;
            }

            cfg.rankList[0][i] = node;
            cfg.rankList[1][i] = random[i];

            cfg.result[i] = node;
            domainResult[i] = cfg.domainList[random[i]];
            numaResult[i] = cfg.numa[node][random[i]];
        }

        // Print the rank file
        printRankFile(nodeList, runTime);
    }

    if (cfg.version == 0)
    {
        // Allocate optimal by using a partition with scotch

        // Step 1:
        // Calculate partitions to reduce inter node communication

        // Calculate amount of processes for each node
        labelList partSizes(cfg.nodes);
        fillAmount(partSizes, cfg.nDomains, cfg.procPN);
        // Info << "partSizes: " << partSizes << endl;

        // if more then one node is used
        if (partSizes.size() > 1)
        {
            decomposeScotch(
                    cfg.optimizeCoeffsDict,
                    "nodePartition",
                    adjncy,
                    xadj,
                    edgeWeights,
                    cfg.procPN,
                    partSizes.size(),
                    partSizes,
                    cfg.result
            );
        }
        else  // all processes on node 0
        {
            cfg.result = 0;
        }

        // Create addressing
        labelList addr(cfg.nDomains);
        labelList counter(cfg.nodes + 1, 0);

        for (label i = 0; i < cfg.nDomains; i++)
        {
            addr[i] = counter[cfg.result[i] + 1];
            counter[cfg.result[i] + 1]++;
        }

        // Step 2:
        // Calculate partitions to reduce inter numa domain communication on each node
        for (label nodei = 0; nodei < cfg.nodes; nodei++)
        {
            labelList subAdjncy, subXadj;
            scalarField subEdgeWeights;

            // build subgraph of node nodei
            subXadj.append(0);
            label count = 0;

            for (label i = 0; i < cfg.result.size(); i++)
            {
                if (cfg.result[i] == nodei)
                {
                    for (label neigh = xadj[i]; neigh < xadj[i+1]; neigh++)
                    {
                        if (cfg.result[adjncy[neigh]] == nodei)
                        {
                            subAdjncy.append(addr[adjncy[neigh]]);
                            subEdgeWeights.append(edgeWeights[neigh]);
                            count++;
                        }
                    }

                    subXadj.append(count);
                }
            }

            label help = cfg.procPN / cfg.numaDomains;

            // Calculate amount of processes for each numa domain
            labelList subPartSizes(cfg.numaDomains, 0);
            fillAmount(subPartSizes, subXadj.size()-1, help);

            labelList subResult(subXadj.size()-1);

            // if more then one numa domain is used
            if (subPartSizes.size() > 1)
            {
                decomposeScotch(
                        cfg.optimizeCoeffsDict,
                        "domainPartition",
                        subAdjncy,
                        subXadj,
                        subEdgeWeights,
                        help,
                        subPartSizes.size(),
                        subPartSizes,
                        subResult
                );
            }
            else  // all processes on numa domain 0
            {
                subResult = 0;
            }

            count = 0;
            for (label i = 0; i < cfg.result.size(); i++)
            {
                if (cfg.result[i] == nodei)
                {
                    domainResult[i] = subResult[count];
                    count++;
                }
            }
        }

        // Create addressing
        labelList subAddr(cfg.nDomains);
        labelListList subCounter(cfg.nodes, labelList(cfg.numaDomains + 1, 0));

        for (label i = 0; i < cfg.nDomains; i++)
        {
            subAddr[i] = subCounter[cfg.result[i]][domainResult[i] + 1];
            subCounter[cfg.result[i]][domainResult[i] + 1]++;
        }

        // Step 3:
        // Calculate partitions to reduce inter numa node communication on each node
        for (label nodei = 0; nodei < cfg.nodes; nodei++)
        {
            for (label domaini = 0; domaini < cfg.numaDomains; domaini++)
            {

                labelList subAdjncy, subXadj;
                scalarField subEdgeWeights;

                // build subgraph of node nodei, numa domain domaini
                subXadj.append(0);
                label count = 0;

                for (label i = 0; i < cfg.nDomains; i++)
                {
                    if (cfg.result[i] == nodei && domainResult[i] == domaini)
                    {
                        for (label neigh = xadj[i]; neigh < xadj[i+1]; neigh++)
                        {
                            if (cfg.result[adjncy[neigh]] == nodei && domainResult[adjncy[neigh]] == domaini)
                            {
                                subAdjncy.append(subAddr[adjncy[neigh]]);
                                subEdgeWeights.append(edgeWeights[neigh]);
                                count++;
                            }
                        }

                        subXadj.append(count);
                    }
                }

                label domainSize = cfg.procPN / (cfg.numaDomains * cfg.numaNodes);

                // Calculate amount of processes for each node
                labelList subPartSizes(cfg.numaNodes * cfg.numaDomains, 0);
                fillAmount(subPartSizes, subXadj.size()-1, domainSize);

                labelList subResult(subXadj.size()-1);

                // if more then one numa node is used
                if (subPartSizes.size() > 1)
                {
                    decomposeScotch(
                            cfg.optimizeCoeffsDict,
                            "numaPartition",
                            subAdjncy,
                            subXadj,
                            subEdgeWeights,
                            cfg.numaNodeSizes,
                            subPartSizes.size(),
                            subPartSizes,
                            subResult
                    );
                }
                else  // all processes on numa node 0
                {
                    subResult = 0;
                }

                count = 0;
                for (label i = 0; i < cfg.result.size(); i++)
                {
                    if (cfg.result[i] == nodei && domainResult[i] == domaini)
                    {
                        numaResult[i] = subResult[count];
                        count++;
                    }
                }
            }
        }

        // sort the new rankfile to start with indizes 0 and sort them ascending
        sortResults(cfg.result, domainResult, numaResult);

        createRankList(cfg.result, domainResult, numaResult);

        // Print the rank file
        printRankFile(nodeList, runTime);
    }

    if (cfg.version == 2)
    {
        // read in cluster process allocation from reportBindings.txt file
        labelList rankToCore(cfg.procPN);
        readBindingFile("reportBindings.txt", rankToCore);

        // create node reults
        for (label i = 0; i < cfg.nodes; i++)
        {
            for (label j = 0; j < cfg.procPN; j++)
            {
                cfg.result[i*cfg.procPN+j] = i;
                cfg.rankList[0][i*cfg.procPN+j] = j;
            }
        }

        // create numa domain and numa nodes results
        for (label i = 0; i < cfg.nodes; i++)
        {
            for (label j = 0; j < cfg.procPN; j++)
            {
                domainResult[i*cfg.procPN+j] = cfg.domainList[rankToCore[j]];
                numaResult[i*cfg.procPN+j] = cfg.numaList[rankToCore[j]];
                cfg.rankList[1][i*cfg.procPN+j] = rankToCore[j];
            }
        }
    }

    calcStatistics(domainResult, numaResult, adjncy, xadj, edgeWeights);

    calcCommunicator(cfg.result, domainResult, numaResult, runTime);

    Info << "\nEnd" << endl;

    }

    return 0;
}




// ************************************************************************* //
