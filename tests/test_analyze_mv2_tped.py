#!/usr/bin/env python

import sys
# For debug, preference local install over all else
if "DEBUG" in sys.argv:
    sys.path.insert(0, "../")
    sys.path.insert(0, "../../") # For mvtest
    sys.path.insert(0, ".")
    sys.argv.remove("DEBUG")

from libgwas.transposed_pedigree_parser import Parser as TransposedPedigreeParser
from libgwas.boundary import BoundaryCheck
from libgwas.snp_boundary_check import SnpBoundaryCheck
from libgwas.data_parser import DataParser
from libgwas.pheno_covar import PhenoCovar
import unittest

import libgwas.standardizer
from meanvar.mvstandardizer import Standardizer
from meanvar import mv_esteq
from . import test_analyze_tped

# We aren't testing the actual application. Just the analysis portion
class TestMv2AnalysisTPed(test_analyze_tped.TestBase):
    def setUp(self):
        super(TestMv2AnalysisTPed, self).setUp()

    def tearDown(self):
        super(TestMv2AnalysisTPed, self).tearDown()

    def testTpedAnalysis(self):
        # We'll start with the correct phenotype with the genotypes, so we'll use
        # a boundary to restrict us to only use the first SNP
        BoundaryCheck.chrom = 1
        DataParser.boundary = BoundaryCheck()
        pheno = PhenoCovar()
        dataset = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        dataset.load_tfam(pheno)
        dataset.load_genotypes()

        results = [x for x in mv_esteq.RunAnalysis(dataset, pheno)]
        self.assertAlmostEqual(0.0034756155, results[0].p_mvtest, places=6)
        self.assertAlmostEqual(0.1134684009, results[0].betas[1], places=6)
        self.assertAlmostEqual(0.0337649965541, results[0].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.0007779211, results[0].beta_pvalues[1], places=6)
        self.assertAlmostEqual(-0.0033479839, results[0].betas[3], places=6)
        self.assertAlmostEqual(0.0492050029324, results[0].beta_stderr[3], places=6)
        self.assertAlmostEqual(0.9457525716, results[0].beta_pvalues[3], places=6)

        self.assertAlmostEqual(0.57778118, results[1].p_mvtest, places=6)
        self.assertAlmostEqual(0.02798537, results[1].betas[1], places=6)
        self.assertAlmostEqual(0.033790691857, results[1].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.40755865, results[1].beta_pvalues[1], places=6)
        self.assertAlmostEqual(0.03275892, results[1].betas[3], places=6)
        self.assertAlmostEqual(0.0475661, results[1].beta_stderr[3], places=6)
        self.assertAlmostEqual(0.49101013, results[1].beta_pvalues[3], places=6)

        self.assertAlmostEqual(0.44661276, results[2].p_mvtest, places=6)
        self.assertAlmostEqual(0.01663975, results[2].betas[1], places=6)
        self.assertAlmostEqual(0.03443300, results[2].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.62891811, results[2].beta_pvalues[1], places=6)
        self.assertAlmostEqual(0.05712017, results[2].betas[3], places=6)
        self.assertAlmostEqual(0.04783608, results[2].beta_stderr[3], places=6)
        self.assertAlmostEqual(0.232446188, results[2].beta_pvalues[3], places=6)


    def testTpedBounded(self):
        BoundaryCheck.chrom = 1
        DataParser.boundary = BoundaryCheck(bp=[2000,3000])
        pheno = PhenoCovar()
        dataset = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        dataset.load_tfam(pheno)
        dataset.load_genotypes()

        results = [x for x in mv_esteq.RunAnalysis(dataset, pheno)]

        self.assertEqual(1, results[0].chr)
        self.assertEqual(2000, results[0].pos)
        self.assertAlmostEqual(0.57778118, results[0].p_mvtest, places=6)
        self.assertAlmostEqual(0.02798537, results[0].betas[1], places=6)
        self.assertAlmostEqual(0.033790691857, results[0].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.40755865, results[0].beta_pvalues[1], places=6)
        self.assertAlmostEqual(0.03275892, results[0].betas[3], places=6)
        self.assertAlmostEqual(0.0475661, results[0].beta_stderr[3], places=6)
        self.assertAlmostEqual(0.49101013, results[0].beta_pvalues[3], places=6)

        self.assertAlmostEqual(0.44661276, results[1].p_mvtest, places=6)
        self.assertAlmostEqual(0.01663975, results[1].betas[1], places=6)
        self.assertAlmostEqual(0.03443300, results[1].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.62891811, results[1].beta_pvalues[1], places=6)
        self.assertAlmostEqual(0.05712017, results[1].betas[3], places=6)
        self.assertAlmostEqual(0.04783608, results[1].beta_stderr[3], places=6)
        self.assertAlmostEqual(0.232446188, results[1].beta_pvalues[3], places=6)


    def testTpedSnpBounded(self):
        BoundaryCheck.chrom = 1
        DataParser.boundary = SnpBoundaryCheck(snps=["rs1000-rs3000"])
        pheno = PhenoCovar()
        dataset = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        dataset.load_tfam(pheno)
        dataset.load_genotypes()


        results = [x for x in mv_esteq.RunAnalysis(dataset, pheno)]
        self.assertEqual(1, results[0].chr)
        self.assertEqual(1000, results[0].pos)

        self.assertAlmostEqual(0.0034756155, results[0].p_mvtest, places=6)
        self.assertAlmostEqual(0.1134684009, results[0].betas[1], places=6)
        self.assertAlmostEqual(0.0337649965541, results[0].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.0007779211, results[0].beta_pvalues[1], places=6)
        self.assertAlmostEqual(-0.0033479839, results[0].betas[3], places=6)
        self.assertAlmostEqual(0.0492050029324, results[0].beta_stderr[3], places=6)
        self.assertAlmostEqual(0.9457525716, results[0].beta_pvalues[3], places=6)

        self.assertAlmostEqual(0.57778118, results[1].p_mvtest, places=6)
        self.assertAlmostEqual(0.02798537, results[1].betas[1], places=6)
        self.assertAlmostEqual(0.033790691857, results[1].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.40755865, results[1].beta_pvalues[1], places=6)
        self.assertAlmostEqual(0.03275892, results[1].betas[3], places=6)
        self.assertAlmostEqual(0.0475661, results[1].beta_stderr[3], places=6)
        self.assertAlmostEqual(0.49101013, results[1].beta_pvalues[3], places=6)

        self.assertAlmostEqual(0.44661276, results[2].p_mvtest, places=6)
        self.assertAlmostEqual(0.01663975, results[2].betas[1], places=6)
        self.assertAlmostEqual(0.03443300, results[2].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.62891811, results[2].beta_pvalues[1], places=6)
        self.assertAlmostEqual(0.05712017, results[2].betas[3], places=6)
        self.assertAlmostEqual(0.04783608, results[2].beta_stderr[3], places=6)
        self.assertAlmostEqual(0.232446188, results[2].beta_pvalues[3], places=6)

class TestMv2AnalysisTPedWithCovariates(test_analyze_tped.TestBase):
    def setUp(self):
        super(TestMv2AnalysisTPedWithCovariates, self).setUp()
        self.standardizer = libgwas.standardizer.get_standardizer()
        libgwas.standardizer.set_standardizer(libgwas.standardizer.NoStandardization)

    def tearDown(self):
        super(TestMv2AnalysisTPedWithCovariates, self).tearDown()
        libgwas.standardizer.set_standardizer(self.standardizer)

    def testTPedAnalysisCov(self):
        PhenoCovar.sex_as_covariate = True
        DataParser.boundary = BoundaryCheck()
        pheno = PhenoCovar()
        dataset = TransposedPedigreeParser(self.tfam_filename, self.tped_filename)
        dataset.load_tfam(pheno)
        dataset.load_genotypes()
        #pheno.standardize_variables()

        results = [x for x in mv_esteq.RunAnalysis(dataset, pheno)]

        self.assertAlmostEqual(0.00342380, results[0].p_mvtest, places=6)
        self.assertAlmostEqual(0.11362883, results[0].betas[1], places=6)
        self.assertAlmostEqual(0.0337610, results[0].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.00076356, results[0].beta_pvalues[1], places=6)
        self.assertAlmostEqual(0.01911090, results[0].betas[3], places=6)
        self.assertAlmostEqual(0.10143178, results[0].beta_stderr[3], places=6)
        self.assertAlmostEqual(0.8505542, results[0].beta_pvalues[3], places=6)


        self.assertAlmostEqual(0.584950593047, results[1].p_mvtest, places=6)
        self.assertAlmostEqual(0.0276543736525, results[1].betas[1], places=6)
        self.assertAlmostEqual(0.03383588, results[1].beta_stderr[1], places=6)
        self.assertAlmostEqual(0.413751829881, results[1].beta_pvalues[1], places=6)

if __name__ == "__main__":
    unittest.main()
