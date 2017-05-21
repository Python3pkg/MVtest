#!/usr/bin/env python
import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../")
    sys.path.insert(0, "../../") # For mvtest
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

import numpy
import os
import mvtest
import unittest
from . import test_analyze_tped
from . import test_analyze_ped
from . import test_pedigree_parser as test_pedigree_parser
from . import test_transped_parser as test_transped_parser
from meanvar import mv_esteq
from libgwas.boundary import BoundaryCheck
from libgwas.data_parser  import DataParser


class TestCmdlineTPed(test_analyze_tped.TestBase):
    def testTPedCmdLineFilenames(self):
        cmds = "--tfile %s" % (self.tfam_filename.split(".")[0])
        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual("PhenoCovar", vars.__class__.__name__)

        self.assertEqual(self.tfam_filename, dataset.tfam_file)
        self.assertEqual(self.tped_filename, dataset.tped_file)
        self.assertEqual(2000, len(dataset.families))
        self.assertEqual([0]*2000, list(dataset.ind_mask[:, 0]))

    def testTPedCmdLineFilenamesExcplicit(self):
        cmds = "--tped %s --tfam %s" % (self.tped_filename, self.tfam_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual("PhenoCovar", vars.__class__.__name__)

        self.assertEqual(self.tfam_filename, dataset.tfam_file)
        self.assertEqual(self.tped_filename, dataset.tped_file)
        self.assertEqual(2000, len(dataset.families))
        self.assertEqual([0]*2000, list(dataset.ind_mask[:, 0]))


    def testTPedCmdLineWithRemove(self):
        missing = "rs3000,rs4000,rs5000"
        cmds = "--tfile %s --exclude %s" % (self.tfam_filename.split(".")[0], missing)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual("Parser", dataset.__class__.__name__)
        self.assertEqual("PhenoCovar", vars.__class__.__name__)

        self.assertEqual(self.tfam_filename, dataset.tfam_file)
        self.assertEqual(self.tped_filename, dataset.tped_file)
        self.assertEqual(2000, len(dataset.families))
        self.assertEqual(missing.split(","), list(DataParser.boundary.ignored_rs))

        results = [x for x in mv_esteq.RunAnalysis(dataset, vars)]
        self.assertEqual(6, len(results))

    def testTPedCmdLineWithExcludeFile(self):
        file = open("__exclusions", "w")
        missing = ["%s:%s" % (i, i) for i in range(0, 500)]
        file.write("\n".join(["%s %s" % (i, i) for i in range(0, 500)]))
        file.close()
        cmds = "--tfile %s --remove __exclusions" % (self.tfam_filename.split(".")[0])

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual("PhenoCovar", vars.__class__.__name__)

        self.assertEqual(self.tfam_filename, dataset.tfam_file)
        self.assertEqual(self.tped_filename, dataset.tped_file)
        self.assertEqual(1500, len(dataset.families))
        self.assertEqual([1]*500 + [0]*1500, list(dataset.ind_mask[:, 0]))
        self.assertEqual(missing, DataParser.ind_exclusions)
        os.remove("__exclusions")

    def testTPedCmdLineWithExclude(self):
        missing = ["%s:%s" % (i, i) for i in range(0, 500)]
        cmds = "--tfile %s --remove %s" % (self.tfam_filename.split(".")[0], ",".join(missing))

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual("PhenoCovar", vars.__class__.__name__)

        self.assertEqual(self.tfam_filename, dataset.tfam_file)
        self.assertEqual(self.tped_filename, dataset.tped_file)
        self.assertEqual(1500, len(dataset.families))
        self.assertEqual([1]*500 + [0]*1500, list(dataset.ind_mask[:, 0]))
        self.assertEqual(",".join(missing), ",".join(DataParser.ind_exclusions))

    def testTPedCmdLineWithBP(self):
        cmds = "--tfile %s --chr=1 --from-bp=1000 --to-bp=5000" % (self.tfam_filename.split(".")[0])

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual(BoundaryCheck.chrom, 1)
        results = [x for x in mv_esteq.RunAnalysis(dataset, vars)]

        self.assertEqual(5, len(results))
        self.assertEqual(1000, results[0].pos)
        self.assertEqual(2000, results[1].pos)
        self.assertEqual(3000, results[2].pos)
        self.assertEqual(4000, results[3].pos)
        self.assertEqual(5000, results[4].pos)

    def testTPedCmdLineWithKB(self):
        cmds = "--tfile %s --chr=1 --from-kb=1 --to-kb=6" % (self.tfam_filename.split(".")[0])

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual(BoundaryCheck.chrom, 1)
        results = [x for x in mv_esteq.RunAnalysis(dataset, vars)]
        self.assertEqual(6, len(results))
        self.assertEqual(1000, results[0].pos)
        self.assertEqual(2000, results[1].pos)
        self.assertEqual(3000, results[2].pos)
        self.assertEqual(4000, results[3].pos)
        self.assertEqual(5000, results[4].pos)
        self.assertEqual(6000, results[5].pos)

    def testTPedCmdLineWithMB(self):
        cmds = "--tfile %s --chr=1 --from-mb=1 --to-mb=2" % (self.tfam_filename.split(".")[0])

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual(1000000, DataParser.boundary.bounds[0])
        self.assertEqual(2000000, DataParser.boundary.bounds[1])

    def testTPedCmdLineWithSNPs(self):
        cmds = "--tfile %s --chr=1 --snps=rs2000-rs4000" % (self.tfam_filename.split(".")[0])

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual(BoundaryCheck.chrom, 1)
        results = [x for x in mv_esteq.RunAnalysis(dataset, vars)]
        self.assertEqual(3, len(results))
        self.assertEqual(2000, results[0].pos)
        self.assertEqual(3000, results[1].pos)
        self.assertEqual(4000, results[2].pos)

    def testTPedCmdLineMAF(self):
        cmds = "--tfile %s --maf=0.3" % (self.tfam_filename.split(".")[0])

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        maf = [0.30225, 0.3075, 0.31, 0.3025, 0.30625]
        i=0
        for snp in dataset:
            self.assertAlmostEqual(maf[i], snp.maf)
            i += 1
    def testTPedCmdLineMaxMAF(self):
        cmds = "--tfile %s --max-maf=0.3" % (self.tfam_filename.split(".")[0])

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        maf = [0.29925, 0.28775, 0.295, 0.2975]
        i=0
        for snp in dataset:
            self.assertAlmostEqual(maf[i], snp.maf)
            i += 1


class TestCmdlinePed(test_analyze_ped.TestBase):
    def testPedCmdLineFilenames(self):
        cmds = "--file %s" % (self.ped_filename.split(".")[0])
        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        self.assertEqual("PhenoCovar", vars.__class__.__name__)

        self.assertEqual(self.map_filename, dataset.mapfile)
        self.assertEqual(self.ped_filename, dataset.datasource)
        self.assertEqual(9, len(dataset.markers))

    def testPedCmdLineMAF(self):
        cmds = "--file %s --maf=0.3" % (self.ped_filename.split(".")[0])

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        maf = [0.30225, 0.3075, 0.31, 0.3025, 0.30625]
        i=0
        for snp in dataset:
            self.assertAlmostEqual(maf[i], snp.maf, places=4)
            i += 1
    def testPedCmdLineMaxMAF(self):
        cmds = "--file %s --max-maf=0.3" % (self.ped_filename.split(".")[0])

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))
        maf = [0.29925, 0.28775, 0.295, 0.2975]
        i=0
        for snp in dataset:
            self.assertAlmostEqual(maf[i], snp.maf)
            i += 1


class TestCmdLineSimplePed(test_pedigree_parser.TestBase):
    def testPedCmdLineMIND(self):
        cmds = "--ped %s --map %s --mind=0.5" % (self.ped_filename_missing, self.map_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        genotypes = [
            [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
            [-1, -1, -1, -1, -1, 1, -1, -1, -1, 0, 0, 1],
            [-1, 2, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [-1, 1, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [-1, 2, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [-1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
             [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]

        mapdata = [x.strip().split() for x in open(self.map_filename).readlines()]

        index = 0
        for snp in dataset:
            self.assertEqual(genotypes[index][1:], list(snp.genotype_data))
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            index += 1
        self.assertEqual(7, index)

    def testPedCmdLineMIND2(self):
        cmds = "--ped %s --map %s --mind=0.10" % (self.ped_filename_missing, self.map_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        genotypes = [
            [ 0, 0, 1, 0],
            [ 1, 0, 0, 1],
            [ 0, 1, 0, 0],
            [ 0, 1, 1, 0],
            [ 0, 1, 0, 0],
            [ 0, 0, 0, 0],
             [0, 0, 0, 0]
        ]

        mapdata = [x.strip().split() for x in open(self.map_filename).readlines()]

        index = 0
        for snp in dataset:
            self.assertEqual(genotypes[index], list(snp.genotype_data))
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            index += 1
        self.assertEqual(5, index)          # Last two are fixed


    def testPedCmdLineGENO(self):
        cmds = "--ped %s --map %s --geno=0.50" % (self.ped_filename_missing, self.map_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        genotypes = [
            [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
            [-1, 2, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [-1, 1, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [-1, 2, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [-1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
             [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]

        mapdata = [['1', 'rs0001', '0', '500'], ['1', 'rs0003', '0', '25000'],
                   ['1', 'rs0004', '0', '45000'], ['2', 'rs0005', '0', '750'],
                   ['2', 'rs0006', '0', '10000'], ['2', 'rs0007', '0', '25000']]

        index = 0
        for snp in dataset:
            self.assertEqual(genotypes[index], list(snp.genotype_data))
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            index += 1
        self.assertEqual(6, index)

    def testPedCmdLineGENO2(self):
        cmds = "--ped %s --map %s --geno=0.05" % (self.ped_filename_missing, self.map_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        genotypes = [
            [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
             [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]

        mapdata = [['1', 'rs0001', '0', '500'], ['2', 'rs0007', '0', '25000']]
        index = 0
        for snp in dataset:
            self.assertEqual(genotypes[index], list(snp.genotype_data))
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            index += 1
        self.assertEqual(2, index)
    def testCmdlineMap3File(self):
        cmds = "--ped %s --map %s --map3 --geno=0.05" % (self.ped_filename_missing, self.map3_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        genotypes = [
            [0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
             [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]

        mapdata = [['1', 'rs0001', '0', '500'], ['2', 'rs0007', '0', '25000']]
        index = 0
        for snp in dataset:
            self.assertEqual(genotypes[index], list(snp.genotype_data))
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            index += 1
        self.assertEqual(2, index)


class TestCmdLineSimpleTPed(test_transped_parser.TestBase):
    def testTPedCmdLineMIND(self):
        cmds = "--tped %s --tfam %s --mind=0.5" % (self.miss_tped_filename, self.tfam_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        genotypes = [
            [0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1],
            [1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]

        mapdata = [['1', 'rs0001', '0', '500'], ['1', 'rs0002', '0', '10000'], ['1', 'rs0003', '0', '25000'],
                   ['1', 'rs0004', '0', '45000'], ['2', 'rs0005', '0', '750'],
                   ['2', 'rs0006', '0', '10000'], ['2', 'rs0007', '0', '25000']]

        index = 0
        for snp in dataset:
            self.assertEqual(genotypes[index], list(snp.genotype_data))
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            index += 1
        self.assertEqual(7, index)

    def testTPedCmdLineMIND2(self):
        cmds = "--tped %s --tfam %s --mind=0.1" % (self.miss_tped_filename, self.tfam_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        genotypes = [
            [0, 1],
            [1, 1],
            [0, 0],
            [0, 0],
            [1, 0],
            [1, 0],
            [0, 0]
        ]

        mapdata = [['1', 'rs0001', '0', '500'], ['1', 'rs0002', '0', '10000'], ['1', 'rs0003', '0', '25000'],
                   ['1', 'rs0004', '0', '45000'], ['2', 'rs0005', '0', '750'],
                   ['2', 'rs0006', '0', '10000'], ['2', 'rs0007', '0', '25000']]

        index = 0
        for snp in dataset:
            self.assertEqual(genotypes[index], list(snp.genotype_data))
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            index += 1
        self.assertEqual(7, index)          # Last two are fixed


    def testTPedCmdLineGENO(self):
        cmds = "--tped %s --tfam %s --geno=0.5" % (self.miss_tped_filename, self.tfam_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        genotypes = [
            [1,  1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
            [0, -1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 0],
            [0, -1, 2, 1, 1, 0, 0, 1, 2, 1, 1, 0],
            [1, -1, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0],
            [1, -1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, -1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0]
        ]

        mapdata = [['1', 'rs0002', '0', '10000'], ['1', 'rs0003', '0', '25000'],
                   ['1', 'rs0004', '0', '45000'], ['2', 'rs0005', '0', '750'],
                   ['2', 'rs0006', '0', '10000'], ['2', 'rs0007', '0', '25000']]

        index = 0
        for snp in dataset:
            self.assertEqual(genotypes[index], list(snp.genotype_data))
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            index += 1
        self.assertEqual(6, index)

    def testTPedCmdLineGENO2(self):
        cmds = "--tped %s --tfam %s --geno=0.05" % (self.miss_tped_filename, self.tfam_filename)

        app = mvtest.MVTestApplication()
        dataset,vars = app.LoadCmdLine(cmds.split(" "))

        genotypes = [
            [1,  1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1]
        ]

        mapdata = [['1', 'rs0002', '0', '10000']]
        index = 0
        for snp in dataset:
            self.assertEqual(genotypes[index], list(snp.genotype_data))
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            index += 1
        self.assertEqual(1, index)
if __name__ == "__main__":
    unittest.main()
