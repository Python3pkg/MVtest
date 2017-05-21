#!/usr/bin/env python

import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../")
    sys.path.insert(0, "../../") # For mvtest
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

from . import test_pedigree_parser as test_pedigree_parser

from libgwas.data_parser import DataParser
from libgwas.pheno_covar import PhenoCovar

import mvtest
import unittest

class TestCmdlineVariedColumns(test_pedigree_parser.TestBase):
    def setUp(self):
        super(TestCmdlineVariedColumns, self).setUp()
        self.sex_as_covar = PhenoCovar.sex_as_covariate
        self.has_sex = DataParser.has_sex
        self.has_parents = DataParser.has_parents
        self.has_fid = DataParser.has_fid
        self.has_pheno = DataParser.has_pheno
        self.has_liability = DataParser.has_liability

    def tearDown(self):
        super(TestCmdlineVariedColumns, self).tearDown()
        PhenoCovar.sex_as_covariate = self.sex_as_covar
        DataParser.has_sex = self.has_sex
        DataParser.has_parents = self.has_parents
        DataParser.has_fid = self.has_fid
        DataParser.has_pheno = self.has_pheno
        DataParser.has_liability = self.has_liability

    def testCmdLinePedigreeNoSex(self):
        f = open(self.ped_filename, "w")
        f.write("""1 1 0 0 0.1 A A G T A A G G C T G T T T
2 2 0 0 0.4 A C G T G G C G T T G G C T
3 3 0 0 1.0 A A G G A G C C C C G T C T
4 4 0 0 0.5 A A G G A G C G C T G G T T
5 5 0 0 0.9 A C G G A A C G C C G G T T
6 6 0 0 1.0 A A G T A A G G C C G G T T
7 7 0 0 0.1 A A G T A A G G C T G G T T
8 8 0 0 0.4 A C G T G G C G T T G G C T
9 9 0 0 1.0 A A G G A G C C C C G T C T
10 10 0 2 0.5 A A G G A G C G C T G G T T
11 11 0 1 0.9 A C G G A A C G C C G G T T
12 12 0 1 1.0 A A G T A A G G C C G G T T""")
        f.close()
        cmds = "--file %s --no-sex" % (self.ped_filename.split(".")[0])
        app = mvtest.MVTestApplication()
        ped_parser,pc = app.LoadCmdLine(cmds.split(" "))

        mapdata = [x.strip().split() for x in open(self.map_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testCmdLinePedigreeWithLiability(self):
        f = open(self.ped_filename, "w")
        f.write("""1 1 0 0 1 0.1 1 A A G T A A G G C T G T T T
2 2 0 0 1 0.4 1 A C G T G G C G T T G G C T
3 3 0 0 2 1.0 1 A A G G A G C C C C G T C T
4 4 0 0 2 0.5 1 A A G G A G C G C T G G T T
5 5 0 0 1 0.9 1 A C G G A A C G C C G G T T
6 6 0 0 1 1.0 1 A A G T A A G G C C G G T T
7 7 0 0 1 0.1 1 A A G T A A G G C T G G T T
8 8 0 0 1 0.4 1 A C G T G G C G T T G G C T
9 9 0 0 2 1.0 1 A A G G A G C C C C G T C T
10 10 0 0 2 0.5 1 A A G G A G C G C T G G T T
11 11 0 0 1 0.9 1 A C G G A A C G C C G G T T
12 12 0 0 1 1.0 1 A A G T A A G G C C G G T T""")

        f.close()
        cmds = "--file %s --liability" % (self.ped_filename.split(".")[0])
        app = mvtest.MVTestApplication()
        ped_parser,pc = app.LoadCmdLine(cmds.split(" "))


        mapdata = [x.strip().split() for x in open(self.map_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testCmdLinePedigreeNoParents(self):
        f = open(self.ped_filename, "w")
        f.write("""1 1 1 0.1 A A G T A A G G C T G T T T
2 2 1 0.4 A C G T G G C G T T G G C T
3 3 2 1.0 A A G G A G C C C C G T C T
4 4 2 0.5 A A G G A G C G C T G G T T
5 5 1 0.9 A C G G A A C G C C G G T T
6 6 1 1.0 A A G T A A G G C C G G T T
7 7 1 0.1 A A G T A A G G C T G G T T
8 8 1 0.4 A C G T G G C G T T G G C T
9 9 2 1.0 A A G G A G C C C C G T C T
10 10 2 0.5 A A G G A G C G C T G G T T
11 11 1 0.9 A C G G A A C G C C G G T T
12 12 1 1.0 A A G T A A G G C C G G T T""")

        f.close()
        cmds = "--file %s --no-parents" % (self.ped_filename.split(".")[0])
        app = mvtest.MVTestApplication()
        ped_parser,pc = app.LoadCmdLine(cmds.split(" "))

        mapdata = [x.strip().split() for x in open(self.map_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testCmdLinePedigreeNoPheno(self):
        f = open(self.ped_filename, "w")
        f.write("""1 1 0 0 1 A A G T A A G G C T G T T T
2 2 0 0 1 A C G T G G C G T T G G C T
3 3 0 0 2 A A G G A G C C C C G T C T
4 4 0 0 2 A A G G A G C G C T G G T T
5 5 0 0 1 A C G G A A C G C C G G T T
6 6 0 0 1 A A G T A A G G C C G G T T
7 7 0 0 1 A A G T A A G G C T G G T T
8 8 0 0 1 A C G T G G C G T T G G C T
9 9 0 0 2 A A G G A G C C C C G T C T
10 10 0 0 2 A A G G A G C G C T G G T T
11 11 0 0 1 A C G G A A C G C C G G T T
12 12 0 0 1 A A G T A A G G C C G G T T""")

        f.close()

        cmds = "--file %s --sex --no-pheno" % (self.ped_filename.split(".")[0])
        app = mvtest.MVTestApplication()
        ped_parser,pc = app.LoadCmdLine(cmds.split(" "))

        mapdata = [x.strip().split() for x in open(self.map_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)
    def testCmdLinePedigreeNoPhenoNoFam(self):
        f = open(self.ped_filename, "w")
        f.write("""1 0 0 1 A A G T A A G G C T G T T T
2 0 0 1 A C G T G G C G T T G G C T
3 0 0 2 A A G G A G C C C C G T C T
4 0 0 2 A A G G A G C G C T G G T T
5 0 0 1 A C G G A A C G C C G G T T
6 0 0 1 A A G T A A G G C C G G T T
7 0 0 1 A A G T A A G G C T G G T T
8 0 0 1 A C G T G G C G T T G G C T
9 0 0 2 A A G G A G C C C C G T C T
10 0 0 2 A A G G A G C G C T G G T T
11 0 0 1 A C G G A A C G C C G G T T
12 0 0 1 A A G T A A G G C C G G T T""")

        f.close()
        cmds = "--file %s --no-fid --no-pheno" % (self.ped_filename.split(".")[0])
        app = mvtest.MVTestApplication()
        ped_parser,pc = app.LoadCmdLine(cmds.split(" "))

        mapdata = [x.strip().split() for x in open(self.map_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)


    def testCmdLinePedigreeNoFamId(self):
        f = open(self.ped_filename, "w")
        f.write("""1 0 0 1 0.1 A A G T A A G G C T G T T T
2 0 0 1 0.4 A C G T G G C G T T G G C T
3 0 0 2 1.0 A A G G A G C C C C G T C T
4 0 0 2 0.5 A A G G A G C G C T G G T T
5 0 0 1 0.9 A C G G A A C G C C G G T T
6 0 0 1 1.0 A A G T A A G G C C G G T T
7 0 0 1 0.1 A A G T A A G G C T G G T T
8 0 0 1 0.4 A C G T G G C G T T G G C T
9 0 0 2 1.0 A A G G A G C C C C G T C T
10 0 0 2 0.5 A A G G A G C G C T G G T T
11 0 0 1 0.9 A C G G A A C G C C G G T T
12 0 0 1 1.0 A A G T A A G G C C G G T T""")

        f.close()
        cmds = "--file %s --no-fid" % (self.ped_filename.split(".")[0])
        app = mvtest.MVTestApplication()
        ped_parser,pc = app.LoadCmdLine(cmds.split(" "))

        mapdata = [x.strip().split() for x in open(self.map_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)

    def testCmdLinePedigreeNoFamSexOrParents(self):
        f = open(self.ped_filename, "w")
        f.write("""1 0.1 A A G T A A G G C T G T T T
2 0.4 A C G T G G C G T T G G C T
3 1.0 A A G G A G C C C C G T C T
4 0.5 A A G G A G C G C T G G T T
5 0.9 A C G G A A C G C C G G T T
6 1.0 A A G T A A G G C C G G T T
7 0.1 A A G T A A G G C T G G T T
8 0.4 A C G T G G C G T T G G C T
9 1.0 A A G G A G C C C C G T C T
10 0.5 A A G G A G C G C T G G T T
11 0.9 A C G G A A C G C C G G T T
12 1.0 A A G T A A G G C C G G T T""")

        f.close()
        cmds = "--file %s --no-parents --no-fid --no-sex" % (self.ped_filename.split(".")[0])
        app = mvtest.MVTestApplication()
        ped_parser,pc = app.LoadCmdLine(cmds.split(" "))

        mapdata = [x.strip().split() for x in open(self.map_filename).readlines()]

        index = 0
        for snp in ped_parser:
            self.assertEqual(int(mapdata[index][0]), snp.chr)
            self.assertEqual(int(mapdata[index][3]), snp.pos)
            self.assertEqual(mapdata[index][1], snp.rsid)
            self.assertEqual(self.genotypes[index], list(snp.genotype_data))
            index += 1
        self.assertEqual(7, index)


if __name__ == "__main__":
    unittest.main()
