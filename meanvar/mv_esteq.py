import sys

import numpy
import scipy
import scipy.stats
import math
import exceptions

from .simple_timer import SimpleTimer
from libgwas.data_parser import DataParser
from .mvresult import MVResult
from libgwas.exceptions import UnsolvedLocus
from libgwas.exceptions import NanInResult
import libgwas.pheno_covar
from libgwas.standardizer import get_standardizer

__copyright__ = "Copyright (C) 2015 Todd Edwards, Chun Li and Eric Torstenson"
__license__ = "GPL3.0"
#     This file is part of MVtest.
#
#     MVtest is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     MVtest is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with MVtest.  If not, see <http://www.gnu.org/licenses/>.

def MeanVarEstEQ(y, x, covariates, tol=1e-8):
    """Perform the mean var calculation using estimated equestions

    :param y: Outcomes
    :param x: [genotypes, cov1, ..., covN]
    :param tol: convergence criterion

    """

    pcount = covariates.shape[0] + 2
    N = y.shape[0]
    beta_count = pcount * 2
    X = numpy.ones((N, pcount))
    X[:, 1] = x

    aprime = [1, 0]
    idx = 2
    for cov in covariates:
        aprime.append(-numpy.mean(cov)/numpy.std(cov))
        X[:, idx] = cov
        idx += 1
    aprime = numpy.matrix(aprime)

    # http://stackoverflow.com/questions/7267226/range-for-floats
    def frange(x, y, jump):
        while x < y:
            yield x
            x += jump
        yield 0

    class PhiReturn(object):
        def __init__(self, phi, dtheta):
            self.phi = phi
            self.dtheta = dtheta

    def dot_diag(x, y):
        """This is a conveniance function to perform the following \
        in a more efficient manner:

        :param x: arr1
        :param y: arr2

        x.dot(numpy.diag(y))

        y must be a single dimensioned array

        """
        if len(y.shape) != 1:
            print("You can't pass an array of shape %s to dot_diag" % (y.shape), file=sys.stderr)
            sys.exit(1)

        result = numpy.empty(x.shape)
        for i in range(0, x.shape[0]):
            result[i] = x[i] * y

        return result

    def diag_dot(x, y):
        """diagonal multiply
        """
        result = numpy.empty(y.shape)
        for i in range(0, y.shape[0]):
            result[i] = y[i] * x[i]
        return result

    def Phi(theta):
        MM = y - numpy.dot(X, theta[0:pcount])
        SS = numpy.exp(numpy.dot(-X, theta[pcount:]))
        DD1 = MM * SS
        DD2 = 0.5 * MM**2 * SS
        Phi = numpy.hstack((DD1.dot(X), (DD2-0.5).dot(X)))

        tX = X.transpose()

        t1 = numpy.empty(tX.shape)
        t2 = numpy.empty(tX.shape)
        t3 = numpy.empty(tX.shape)
        for i in range(0,tX.shape[0]):
            t1[i] = tX[i]*-SS
            t2[i] = tX[i]*-DD1
            t3[i] = tX[i]*-DD2
        t1 = t1.dot(X)
        t2 = t2.dot(X)
        t3 = t3.dot(X)

        dtheta = numpy.hstack((numpy.vstack((t1, t2)), numpy.vstack((t2, t3))))
        return PhiReturn(Phi, dtheta)



    class MvSolveReturn(object):
        def __init__(self, theta, dtheta):
            self.theta = theta
            self.dtheta = dtheta

    def MVsolve(theta_new):
        solution_found = False

        mvsolve_iterations = 0

        while not solution_found:
            theta_old = theta_new.copy()
            tmp = Phi(theta_old)
            #print "ITR", mvsolve_iterations, numpy.mean(theta_new),numpy.sum(numpy.absolute(theta_new - theta_old))
            theta_new = theta_old - scipy.linalg.solve(tmp.dtheta, tmp.phi)
            mvsolve_iterations += 1

            if (numpy.sum(numpy.absolute(theta_new - theta_old)) < tol):
                tmp = Phi(theta_new)
                solution_found = True
            else:
                if mvsolve_iterations > 25000:
                    #print >> sys.stderr, mvsolve_iterations, "failures"
                    raise UnsolvedLocus("")

        return MvSolveReturn(theta_new, tmp.dtheta), mvsolve_iterations

    def MVcalcB(theta):
        MM = y - X.dot(theta[0:pcount])
        SS = numpy.exp(-X.dot(theta[pcount:]))
        DD1 = MM * SS
        DD2 = 0.5 * MM**2 * SS
        AA = numpy.hstack((diag_dot(DD1, X), diag_dot(DD2-0.5,X)))

        return numpy.transpose(AA).dot(AA)/ N


    mod = None
    itr = 0
    total_iterations = 0
    for i in frange(0.00, 1.0, 0.05):
        theta=numpy.empty((beta_count))
        theta[:] = i

        try:
            mod, iterations = MVsolve(theta)
            total_iterations += iterations
            if i > 0.05:
                print("Completed: ", total_iterations, itr, file=sys.stderr)
            break
        except exceptions.ValueError as e:
            pass
        except numpy.linalg.linalg.LinAlgError as e:
            pass
        except Exception as inst:
            #print type(inst)
            pass
        itr += 1

    if not mod:
        raise UnsolvedLocus("")

    try:
        ainv = scipy.linalg.inv(mod.dtheta) * N
    except:
        raise ValueError("Singular Matrix Encountered")
    B = MVcalcB(mod.theta)
    V = ainv.dot(B).dot(ainv.transpose())



    # Focus on the two parameters of interest
    theta2 = numpy.array([mod.theta[1], mod.theta[pcount+1]])
    V2 = V[1:beta_count:pcount,1:beta_count:pcount]
    pvalt = 1 - scipy.stats.chi2.cdf(theta2.dot(scipy.linalg.inv(V2)).dot(theta2) * N, 2)


    ## From Chun's updated code:
    theta= mod.theta
    se = numpy.sqrt(numpy.diag(V)/N)
    pval = 2*scipy.stats.norm.cdf(-numpy.absolute(mod.theta/se))
    return pvalt, theta, pval, se, V/N

def RunMeanVar(pheno, geno, covar=[]):
    """Setup and execute the mean var calculation.

    :param pheno: Phenotype data (one phenotype at a time)
    :param geno: SNP data (might be genotypes, or dosages, etc)
    :param covar: List of covariate data

    It is possible that the optimization will fail to converge. Such
    cases are stripped of data, but are still reported to alert the
    user that there were problems with the data.

    """
    return MeanVarEstEQ(pheno, geno, covar)

def RunAnalysis(dataset, pheno_covar):
    """Run the actual analysis on all valid loci for each phenotype

    :param dataset: GWAS parser object
    :param pheno_covar: holds all of the variables

    This acts as a standard iterator, returning a single MVResult for
    each locus/phenotype combination.

    Missing is evaluated as anything missing in any of the
    phenotype, covariate(s) or genotype

    """
    covar_missing = numpy.empty(dataset.ind_count, dtype=bool)
    covar_missing[:] = False

    total_covar_count = len(pheno_covar.covariate_data)
    iteration_count = []

    unsolved = []

    pcount = 2+total_covar_count

    std = get_standardizer()

    for snp in dataset:
        for y in pheno_covar:
            st = SimpleTimer()
            (pheno, covariates, nonmissing) = y.get_variables((snp.genotype_data==DataParser.missing_storage))


            genotypes  = snp.genotype_data[nonmissing]

            try:
                pvalt, estimates, pvalues, se, v = RunMeanVar(pheno, genotypes, covariates)

                lmgeno = genotypes
                for c in covariates:
                    lmgeno = lmgeno + c
                lm = scipy.stats.linregress(pheno, lmgeno)[3]


                betavars, se, pvalues = y.destandardize(estimates, se, pvalues=pvalues,v=v, nonmissing=nonmissing)

                betas = betavars[0:pcount]
                betase = se[0:pcount]
                vars = betavars[pcount:]
                varse = se[pcount:]

                # Check for invalid output
                for var in betavars + se + list(pvalues):
                    if numpy.isnan(var):
                        raise NanInResult()

                nonmissing_ct = numpy.sum(nonmissing)

                result = MVResult(
                                snp.chr,
                                snp.pos,
                                snp.rsid,
                                snp.major_allele,
                                snp.minor_allele,
                                dataset.get_effa_freq(genotypes),
                                non_miss_count=nonmissing_ct,
                                ph_label=y.get_phenotype_name(),
                                p_mvtest=pvalt,
                                beta_values=list(betas)+list(vars),
                                pvalues=pvalues,
                                stderrors=list(betase) + list(varse),
                                maf=snp.maf,
                                covar_labels=y.get_covariate_names(),
                                lm=lm,
                                runtime = st.runtime(),
                )

                result.blin = estimates[0:pcount]
                result.bvar = estimates[pcount:]
                yield result
            except NanInResult as e:
                print("\t".join([str(x) for x in [
                                snp.chr,
                                snp.pos,
                                snp.rsid,
                                y.get_phenotype_name(),
                                "%d" % (numpy.sum(nonmissing)),
                                snp.major_allele,
                                snp.minor_allele,
                                snp.allele_count2,
                                "NAN-Found",
                                "MAF=%0.4f" % (snp.maf)]]), file=sys.stderr)
            except ValueError as e:
                print("\t".join([str(x) for x in [
                                snp.chr,
                                snp.pos,
                                snp.rsid,
                                y.get_phenotype_name(),
                                "%d" % (numpy.sum(nonmissing)),
                                snp.major_allele,
                                snp.minor_allele,
                                snp.allele_count2,
                                "Unsolvable",
                                "MAF=%0.4f" % (snp.maf)]]), file=sys.stderr)
            except UnsolvedLocus as e:
                print("\t".join([str(x) for x in [
                                snp.chr,
                                snp.pos,
                                snp.rsid,
                                y.get_phenotype_name(),
                                "%d" % (numpy.sum(nonmissing)),
                                snp.major_allele,
                                snp.minor_allele,
                                snp.allele_count2,
                                "Unsolved",
                                "MAF=%0.4f" % (snp.maf)]]), file=sys.stderr)
                unsolved.append(snp)
    if len(unsolved)>0:
        print("Total unsolvable loci: ", len(unsolved), file=sys.stderr)

