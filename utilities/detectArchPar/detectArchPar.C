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
    Sends point-to-point communication messages between each pair of
    used processes in the OpenFOAM case bidirectionally and measures
    the time on the sending side. The measured times are print as a
    matrix in a csv file for further processing.

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
#include "processorFvPatchField.H"

#include <asm-generic/errno.h>
#include <cwctype>
#include <fstream>
#include <cstdio>
#include <pthread.h>
#include <regex>
#include <algorithm>
#include <cstdio>   // For popen()
#include <cstdlib>  // For atoi()


#include <iostream>
#include <chrono>
#include <unistd.h>
#include <sstream>
#include <iomanip>

#include <mpi.h>


extern "C"
{
#include "scotch.h"
}

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct NumaConfig
{
    label rank;
    label nDomains;
    label nodes;
    label procPN;
    label numaNodeSizes;
    label numaDomains;
    label numaNodes;

    label version;

    dictionary decomposeParDict;
    dictionary optimizeCoeffsDict;

    labelList numaList;
    labelList domainList;

    // Default constructor
    NumaConfig()
    :
        rank(0),
        nDomains(0),
        nodes(0),
        procPN(0),
        numaNodeSizes(0),
        numaDomains(0),
        numaNodes(0),
        version(0),
        decomposeParDict(),
        optimizeCoeffsDict(),
        numaList(),
        domainList()
    {}
};

NumaConfig cfg;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

void test(const label& testSize, const label runs, const bool& writeTime)
{
    // Initialize scalar send and receive data according to the testSize
    Field<scalar> sendData(testSize, 123.456);
    Field<scalar> recvData(testSize, 0.0);

    Field<scalar> resetField(testSize, 0.0);

    label comm = 0;

    // Lists for directly calculating the mean times per comm. layer
    scalarList erg(4, 0.0);
    labelList ergCount(4, 0);

    // List for measured times and final mean time matrix
    scalarList time(cfg.nDomains, 0.0);
    scalarListList timeMat(cfg.nDomains, scalarList(cfg.nDomains, 0.0));

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();

    {
        // Execute test for a given amount of test runs
        for (label run = 0; run < runs; run++)
        {
            for (label send = 0; send < cfg.nDomains; send++)
            {
                for (label recv = 0; recv < cfg.nDomains; recv++)
                {
                    if(send == recv) continue;

                    // Let each pair of processes communiacte bidirectionally
                    // Info << "send: " << send << " -> " << recv << " :recv\n";

                    if (cfg.rank == recv)
                    {
                        UIPstream::read
                        (
                            Pstream::commsTypes::blocking,
                            send, // source rank
                            reinterpret_cast<char*>(recvData.begin()),
                            recvData.byteSize(),
                            recv,
                            comm
                        );

                        UOPstream::write
                        (
                            Pstream::commsTypes::blocking,
                            send, // destination rank
                            reinterpret_cast<const char*>(recvData.begin()),
                            recvData.byteSize(),
                            recv,
                            comm
                        );
                    }

                    if (cfg.rank == send)
                    {
                        start = std::chrono::high_resolution_clock::now();

                        UOPstream::write
                        (
                            Pstream::commsTypes::blocking,
                            recv, // destination rank
                            reinterpret_cast<const char*>(sendData.begin()),
                            sendData.byteSize(),
                            recv,
                            comm
                        );

                        UIPstream::read
                        (
                            Pstream::commsTypes::blocking,
                            recv, // source rank
                            reinterpret_cast<char*>(recvData.begin()),
                            recvData.byteSize(),
                            recv,
                            comm
                        );

                        end = std::chrono::high_resolution_clock::now();
                        std::chrono::duration<double, std::micro> duration = end - start;

                        if (recvData[0] != 123.456)
                        {
                            Pout << "Wrong or failed communication between ranks "
                                << send << " and " << recv << "!" << endl;
                        }

                        // Use half the measured time as the one way comm. time
                        time[recv] += duration.count() / 2;

                        recvData = resetField;
                    }
                }
            }
        }

        // Calculate mean of the measured time values
        for (label i = 0; i < time.size(); i++)
        {
            time[i] /= runs;
        }

        // Collect all measured mean times on first proc of the node
        if (cfg.rank == 0)
        {
            for (label i = 0; i < cfg.nDomains; i++)
            {
                timeMat[0][i] = time[i];
            }
        }

        for (label proci = 1; proci < cfg.nDomains; proci++)
        {
            if (cfg.rank == 0)
            {
                label recvReq = UPstream::nRequests();
                UIPstream::read
                (
                    Pstream::commsTypes::nonBlocking,
                    proci, // source rank
                    reinterpret_cast<char*>(time.begin()),
                    time.byteSize(),
                    proci,
                    comm
                );

                if(recvReq >= 0)
                {
                    UPstream::waitRequest(recvReq);
                }
            }

            // send measured time values to rank 0
            if (cfg.rank == proci)
            {
                UOPstream::write
                (
                    Pstream::commsTypes::nonBlocking,
                    0, // destination rank
                    reinterpret_cast<const char*>(time.begin()),
                    time.byteSize(),
                    proci,
                    comm
                );
            }

            for (label i = 0; i < cfg.nDomains; i++)
            {
                timeMat[proci][i] = time[i];
            }
        }

        UPstream::resetRequests(0);
    }

    // Print time matrix to log file
    Info << "Time Matrix:" << endl;
    for (label i = 0; i < timeMat.size(); i++)
    {
        for (label j = 0; j < timeMat[i].size(); j++)
        {
            Info << timeMat[i][j] << "  ";
        }
        Info << endl;
    }

    // Write time matrix into csv file
    if (writeTime)
    {
        if (cfg.rank == 0)
        {
            string filename = "./timeMat.csv";

            std::ofstream out;
            out.open(filename, std::ios::app);

            if (!out) {
                    Foam::FatalError << "Could not open file "
                    << filename << nl
                    << Foam::abort(Foam::FatalError);
            }

            out << std::fixed << std::setprecision(4);

            for (label i = 0; i < timeMat.size(); i++)
            {
                for (label j = 0; j < timeMat[i].size()-1; j++)
                {
                    out << timeMat[i][j] << ",";
                }
                out << timeMat[i][timeMat[i].size()-1];
                out << "\n";
            }

            out.close();
        }
    }

    // Calculate mean times per comm. layer according to given hardware topology
    for (label send = 0; send < cfg.nDomains; send++)
    {
        for (label recv = 0; recv < cfg.nDomains; recv++)
        {
            if(send == recv) continue;

            if(timeMat[send][recv] > 1000) continue;

            if(abs(send-recv) <= 1) continue;

            if (send / cfg.procPN != recv / cfg.procPN)
            {
                erg[3] += timeMat[send][recv];
                ergCount[3]++;
            }
            else
            {
                label sendSlot = send % cfg.procPN;
                label recvSlot = recv % cfg.procPN;

                if (cfg.domainList[sendSlot] != cfg.domainList[recvSlot])
                {
                    erg[2] += timeMat[send][recv];
                    ergCount[2]++;
                }
                else if (cfg.numaList[sendSlot] != cfg.numaList[recvSlot])
                {
                    erg[1] += timeMat[send][recv];
                    ergCount[1]++;
                }
                else
                {
                    erg[0] += timeMat[send][recv];
                    ergCount[0]++;
                }
            }
        }
    }


    for (label i = 0; i < erg.size(); i++)
    {
        if (ergCount[i] != 0)
        {
            erg[i] /= ergCount[i];
        }
    }


    Info << "\nCommunication Test with " << runs << " Testruns and a Messagesize " << testSize << " for\n";
    Info << "(intra-NUMA-node, inter-NUMA-node, inter-NUMA-domain, inter-node)\n";
    Info << erg << endl << endl;

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
            << "detectArchPar must be run with -parallel"
            << exit(FatalError);
    }

    cfg.rank = Pstream::myProcNo();

    // Set time from database
    processorRunTimes runTimes(Foam::Time::controlDictName, args);

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

    // Get size and iterations of test communication message
    label testSize = cfg.optimizeCoeffsDict.lookupOrDefault<label>("testSize", 2000);
    label testIterations = cfg.optimizeCoeffsDict.lookupOrDefault<label>("testIterations", 1000);

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
    Info << "testSize: " << testSize << endl;
    Info << "testIterations: " << testIterations << endl;

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

    // Read in numaList and domainList (for all versions -> calc statistics)
    {
        // Get the numa list from the subdictionary
        cfg.numaList = cfg.optimizeCoeffsDict.lookupOrDefault<labelList>("numaList", labelList());

        Info << "numaList:" << cfg.numaList << endl;

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

        Info << "domainList: " << cfg.domainList << endl;

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

    bool writeTimeMat = true;

    Info << "Test communication ..." << endl;
    test(testSize, testIterations, writeTimeMat);
    Info << "DONE!" << endl << endl;

    Info << "\nEnd" << endl;

    return 0;
}




// ************************************************************************* //
