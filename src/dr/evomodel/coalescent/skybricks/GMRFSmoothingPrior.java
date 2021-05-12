package dr.evomodel.coalescent.skybricks;


import dr.inference.glm.GeneralizedLinearModel;
import dr.inference.model.*;
import dr.xml.XMLParseException;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;

import java.util.ArrayList;
import java.util.List;

public interface GMRFSmoothingPrior extends Model, Likelihood {
    public static final double LOG_TWO_TIMES_PI = 1.837877;

    double getLogFieldLikelihood();

    DenseVector getMeanAdjustedGamma();

    SymmTridiagMatrix getStoredScaledWeightMatrix(double precision, double lambda);

    SymmTridiagMatrix getScaledWeightMatrix(double precision, double lambda);

    void setupGMRFWeights();

    void setupGMRFWeights(double[] gridPoints, double fieldScalar);

    Parameter getLambdaParameter();

    SymmTridiagMatrix getWeightMatrix();

    Parameter getPrecisionParameter();

    Parameter getSmoothedParameter();

    Parameter getBetaParameter();

    MatrixParameter getDesignMatrix();

    List<Parameter> getBetaListParameter();

    List<MatrixParameter> getCovariates();

    double[] getGradientWrtPrecision();

    double[] getDiagonalHessianWrtPrecision();

    double[] getDiagonalHessianWrtRegressionCoefficients();

    double[] getGradientWrtRegressionCoefficients();

    /**
     * This builder class hides the nasty logic about what smoothing prior is needed given the covariates, glm model,
     * or lack there of.
     */
    class Builder {
        //TODO replace fieldLength with epochModel
        private Parameter smoothingParameter;
        private Parameter precisionParameter;
        private Parameter lambdaParameter;
        private List<Parameter> betaList = null;
        private MatrixParameter dMatrix = null;
        private List<MatrixParameter> covariates = null;
        private Parameter betaParameter = null;
        private List<Parameter> deltaList = null;
        private List<Parameter> covPrecParametersRecent = null;
        private List<Parameter> covPrecParametersDistant = null;
        private List<Parameter> firstObservedIndexParameter = null;
        private List<Parameter> lastObservedIndexParameter = null;
        private Parameter recentIndices = null;
        private Parameter distantIndices = null;
        private EpochProvider epochModel = null;
        private boolean useGlmModel = false;
        private GeneralizedLinearModel glm = null;

        public Builder() {
        }


        public void setSmoothingParameter(Parameter smoothingParameter) {
            this.smoothingParameter = smoothingParameter;
        }

        public void setPrecisionParameter(Parameter precisionParameter) {
            this.precisionParameter = precisionParameter;
        }

        public void setLambdaParameter(Parameter lambdaParameter) {
            this.lambdaParameter = lambdaParameter;
        }

        public void setBetaList(List<Parameter> betaList) {
            this.betaList = betaList;
        }

        public void setdMatrix(MatrixParameter dMatrix) {
            this.dMatrix = dMatrix;
        }

        public void setCovariates(List<MatrixParameter> covariates) {
            this.covariates = covariates;
        }

        public void setBetaParameter(Parameter betaParameter) {
            this.betaParameter = betaParameter;
        }

        public void setDeltaList(List<Parameter> deltaList) {
            this.deltaList = deltaList;
        }

        public void setCovPrecParametersRecent(List<Parameter> covPrecParametersRecent) {
            this.covPrecParametersRecent = covPrecParametersRecent;
        }

        public void setCovPrecParametersDistant(List<Parameter> covPrecParametersDistant) {
            this.covPrecParametersDistant = covPrecParametersDistant;
        }

        public void setFirstObservedIndexParameter(List<Parameter> firstObservedIndexParameter) {
            this.firstObservedIndexParameter = firstObservedIndexParameter;
        }

        public void setLastObservedIndexParameter(List<Parameter> lastObservedIndexParameter) {
            this.lastObservedIndexParameter = lastObservedIndexParameter;
        }

        public void setRecentIndices(Parameter recentIndices) {
            this.recentIndices = recentIndices;
        }

        public void setDistantIndices(Parameter distantIndices) {
            this.distantIndices = distantIndices;
        }

        public void setEpochModel(EpochProvider epochModel) {
            this.epochModel = epochModel;
        }

        public void setUseGlmModel(boolean useGlmModel) {
            this.useGlmModel= useGlmModel;
        }

        public void setGlm(GeneralizedLinearModel glm) {
            this.glm = glm;
        }



        public GMRFSmoothingPrior build() throws XMLParseException {

            if (dMatrix != null) {
                if (dMatrix.getRowDimension() != smoothingParameter.getDimension())
                    throw new XMLParseException("Design matrix row dimension must equal the smoothed parameter length.");
                if (dMatrix.getColumnDimension() != betaParameter.getDimension())
                    throw new XMLParseException("Design matrix column dimension must equal the regression coefficient length.");
            }

            if((covPrecParametersDistant == null && lastObservedIndexParameter != null) || (covPrecParametersDistant != null && lastObservedIndexParameter == null)){
                throw new XMLParseException("Must specify both lastObservedIndex and covariatePrecision");
            }

            if((covPrecParametersRecent == null && firstObservedIndexParameter != null) || (covPrecParametersRecent != null && firstObservedIndexParameter == null)){
                throw new XMLParseException("Must specify both firstObservedIndex and covariatePrecision");
            }

            if(firstObservedIndexParameter == null && recentIndices != null){
                throw new XMLParseException("Cannot specify covIndicesMissingRecent without specifying firstObservedIndex");
            }
            if(lastObservedIndexParameter == null && distantIndices != null){
                throw new XMLParseException("Cannot specify covIndicesMissingDistant without specifying lastObservedIndex");
            }

            if ((covariates != null && betaList == null) ||
                    (covariates == null &&  betaList != null))
                throw new XMLParseException("Must specify both a set of regression coefficients and a design matrix.");

            if(useGlmModel) {

                covariates = new ArrayList<MatrixParameter>();
                betaList = new ArrayList<Parameter>();
                List<DesignMatrix> designMat = glm.getDesignMatrix();
                List<Parameter> indepParam = glm.getIndependentParameter();
                List<Parameter> indepParamDelta = glm.getIndependentParameterDelta();
                deltaList = new ArrayList<Parameter>();

                for(int i = 0; i < indepParam.get(0).getSize(); i++){
                    MatrixParameter matParam = new MatrixParameter("covariate values", 1,designMat.get(0).getRowDimension());

                    for(int j = 0; j < matParam.getRowDimension(); j++){
                        matParam.setParameterValue(0, j, designMat.get(0).getParameterValue(0, j));
                    }
                    covariates.add(matParam);

                    Parameter betaParam = new Parameter.Default(1);
                    betaParam.setParameterValue(0, indepParam.get(0).getParameterValue(i));
                    betaList.add(betaParam);

                    if(indepParamDelta != null){
                        Parameter deltaParam = new Parameter.Default(1);
                        deltaParam.setParameterValue(0, indepParamDelta.get(0).getParameterValue(i));
                        deltaList.add(deltaParam);
                    }
                }

            }

            if (betaList != null || betaParameter != null) {
                if (lastObservedIndexParameter != null || firstObservedIndexParameter != null) {
                    return new MissingCovariates( smoothingParameter, precisionParameter, lambdaParameter,
                            betaList, dMatrix, covariates, betaParameter, deltaList, covPrecParametersRecent, covPrecParametersDistant, firstObservedIndexParameter, lastObservedIndexParameter, recentIndices, distantIndices);
                } else {
                    return new Covariates( smoothingParameter, precisionParameter, lambdaParameter,betaList,
                            dMatrix, covariates, betaParameter, deltaList);
                }
            } else {
                return new Default( smoothingParameter, precisionParameter, lambdaParameter);
            }
        }
    }

    class Default extends AbstractModelLikelihood implements GMRFSmoothingPrior {


        public Default( Parameter smoothedParameter, Parameter precisionParameter, Parameter lambdaParameter, Parameter betaParameter, MatrixParameter dMatrix) {
            super("GMRFsmoothingPrior");

            this.fieldLength = smoothedParameter.getDimension();
            this.smoothedParameter = smoothedParameter;
            this.precisionParameter = precisionParameter;
            this.lambdaParameter = lambdaParameter;
            this.dMatrix = dMatrix;
            this.betaParameter = betaParameter;

            this.likelihoodKnown = false;

            addVariable(smoothedParameter);
            addVariable(precisionParameter);
            addVariable(lambdaParameter);

            if (betaParameter != null) {
                addVariable(betaParameter);
            }
            if (dMatrix != null) {
                addVariable(dMatrix);
            }
            setupGMRFWeights();
        }

        public Default( Parameter smoothedParameter, Parameter precisionParameter, Parameter lambdaParameter) {
            this(smoothedParameter, precisionParameter, lambdaParameter, null, null);
        }

        public SymmTridiagMatrix getWeightMatrix() {
            return weightMatrix.copy();
        }

        public Parameter getSmoothedParameter() {
            return this.smoothedParameter;
        }

        ;

        public Parameter getLambdaParameter() {
            return this.lambdaParameter;
        }

        public Parameter getPrecisionParameter() {
            return this.precisionParameter;
        }

        public Parameter getBetaParameter() {
            return this.betaParameter;
        }

        public MatrixParameter getDesignMatrix() {
            return this.dMatrix;
        }

        public List<Parameter> getBetaListParameter() {
            return null;
        }

        @Override
        public List<MatrixParameter> getCovariates() {
            return null;
        }

        public DenseVector getMeanAdjustedGamma() {
            DenseVector currentGamma = new DenseVector(smoothedParameter.getAttributeValue());
            updateGammaWithCovariates(currentGamma);
            return currentGamma;
        }

        double handleMissingValues() {
            return 0.0;
        }

        void updateGammaWithCovariates(DenseVector currentGamma) {
            // Do nothing
        }

        public void setupGMRFWeights(double[] gridPoints, double fieldScalar) {
            //Set up the weight Matrix
            double[] offdiag = new double[fieldLength - 1];
            double[] diag = new double[fieldLength];

            //First set up the offdiagonal entries;

            for (int i = 0; i < fieldLength - 1; i++) {
                offdiag[i] = -2.0 / (gridPoints[i] + gridPoints[i + 1]) * fieldScalar;
            }


            //Then set up the diagonal entries;
            for (int i = 1; i < fieldLength - 1; i++)
                diag[i] = -(offdiag[i] + offdiag[i - 1]);

            //Take care of the endpoints
            diag[0] = -offdiag[0];
            diag[fieldLength - 1] = -offdiag[fieldLength - 2];

            weightMatrix = new SymmTridiagMatrix(diag, offdiag);

        }

        public void setupGMRFWeights() {

            //Set up the weight Matrix
            double[] offdiag = new double[fieldLength - 1];
            double[] diag = new double[fieldLength];

            double diagonalValue = 2;
            //First set up the offdiagonal entries;

            for (int i = 0; i < fieldLength - 1; i++) {
                offdiag[i] = -1;
            }

            //Then set up the diagonal entries;
            for (int i = 1; i < fieldLength - 1; i++) {
                //	diag[i] = -(offdiag[i] + offdiag[i - 1]);
                diag[i] = diagonalValue;
            }
            //Take care of the endpoints
            //diag[0] = -offdiag[0];
            //diag[fieldLength - 1] = -offdiag[fieldLength - 2];
            diag[0] = diagonalValue - 1.0;
            diag[fieldLength - 1] = diagonalValue - 1.0;

            weightMatrix = new SymmTridiagMatrix(diag, offdiag);

        }


        public SymmTridiagMatrix getScaledWeightMatrix(double precision) {
            SymmTridiagMatrix a = weightMatrix.copy();
            for (int i = 0; i < a.numRows() - 1; i++) {
                a.set(i, i, a.get(i, i) * precision);
                a.set(i + 1, i, a.get(i + 1, i) * precision);
            }
            a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
            return a;
        }


        public SymmTridiagMatrix getScaledWeightMatrix(double precision, double lambda) {
            if (lambda == 1)
                return getScaledWeightMatrix(precision);

            SymmTridiagMatrix a = weightMatrix.copy();
            for (int i = 0; i < a.numRows() - 1; i++) {
                a.set(i, i, precision * (1 - lambda + lambda * a.get(i, i)));
                a.set(i + 1, i, a.get(i + 1, i) * precision * lambda);
            }

            a.set(fieldLength - 1, fieldLength - 1, precision * (1 - lambda + lambda * a.get(fieldLength - 1, fieldLength - 1)));
            return a;
        }

        public SymmTridiagMatrix getStoredScaledWeightMatrix(double precision) {
            SymmTridiagMatrix a = weightMatrix.copy();
            for (int i = 0; i < a.numRows() - 1; i++) {
                a.set(i, i, a.get(i, i) * precision);
                a.set(i + 1, i, a.get(i + 1, i) * precision);
            }
            a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
            return a;
        }

        public SymmTridiagMatrix getStoredScaledWeightMatrix(double precision, double lambda) {
            if (lambda == 1)
                return getStoredScaledWeightMatrix(precision);

            SymmTridiagMatrix a = storedWeightMatrix.copy();
            for (int i = 0; i < a.numRows() - 1; i++) {
                a.set(i, i, precision * (1 - lambda + lambda * a.get(i, i)));
                a.set(i + 1, i, a.get(i + 1, i) * precision * lambda);
            }

            a.set(fieldLength - 1, fieldLength - 1, precision * (1 - lambda + lambda * a.get(fieldLength - 1, fieldLength - 1)));
            return a;
        }


        public double[] getGradientWrtPrecision() {
            double[] gamma = getMeanAdjustedGamma().getData();
            double currentPrec = getPrecisionParameter().getParameterValue(0);
            int numGridPoints = fieldLength - 1;
            double grad = numGridPoints / (2 * currentPrec);
            for (int i = 0; i < numGridPoints; i++) {
                grad -= 0.5 * (gamma[i + 1] - gamma[i]) * (gamma[i + 1] - gamma[i]);
            }

            return new double[]{grad};
        }

        @Override
        public double[] getDiagonalHessianWrtPrecision() {
            int numGridPoints = fieldLength - 1;

            double currentPrec = getPrecisionParameter().getParameterValue(0);
            double hessian = -numGridPoints / (2 * currentPrec * currentPrec);

            return new double[]{hessian};
        }

        @Override
        public double[] getDiagonalHessianWrtRegressionCoefficients() {
            return null;
        }

        @Override
        public double[] getGradientWrtRegressionCoefficients() {
            return null;
        }


        public double getLogFieldLikelihood() {

            DenseVector diagonal1 = new DenseVector(fieldLength);
            DenseVector currentGamma = getMeanAdjustedGamma();

            double currentLike = handleMissingValues();

            SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
            currentQ.mult(currentGamma, diagonal1);

            currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
            if (lambdaParameter.getParameterValue(0) == 1) {
                currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
            } else {
                currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
            }

            return currentLike;
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            likelihoodKnown = false;
        }

        /**
         * This method is called whenever a parameter is changed.
         * <p/>
         * It is strongly recommended that the model component sets a "dirty" flag and does no
         * further calculations. Recalculation is typically done when the model component is asked for
         * some information that requires them. This mechanism is 'lazy' so that this method
         * can be safely called multiple times with minimal computational cost.
         *
         * @param variable
         * @param index
         * @param type
         */
        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            likelihoodKnown = false;
        }

        /**
         * Additional state information, outside of the sub-model is stored by this call.
         */
        @Override
        protected void storeState() {

            storedWeightMatrix = weightMatrix.copy();
            storedLikelihoodKnown = likelihoodKnown;
            storedLogLikelihood = logLikelihood;

        }

        /**
         * After this call the model is guaranteed to have returned its extra state information to
         * the values coinciding with the last storeState call.
         * Sub-models are handled automatically and do not need to be considered in this method.
         */
        @Override
        protected void restoreState() {

            weightMatrix = storedWeightMatrix;
            likelihoodKnown = storedLikelihoodKnown;
            logLikelihood = storedLogLikelihood;
        }

        /**
         * This call specifies that the current state is accept. Most models will not need to do anything.
         * Sub-models are handled automatically and do not need to be considered in this method.
         */
        @Override
        protected void acceptState() {

        }

        /**
         * Get the model.
         *
         * @return the model.
         */
        @Override
        public Model getModel() {
            return this;
        }

        /**
         * Get the log likelihood.
         *
         * @return the log likelihood.
         */
        @Override
        public double getLogLikelihood() {
            if (!likelihoodKnown) {
                logLikelihood = getLogFieldLikelihood();
            }
            return logLikelihood;
        }

        /**
         * Forces a complete recalculation of the likelihood next time getLikelihood is called
         */
        @Override
        public void makeDirty() {
            likelihoodKnown = false;
        }

        private double storedLogLikelihood;
        private boolean likelihoodKnown;
        private boolean storedLikelihoodKnown;
        private double logLikelihood;

        int fieldLength;
        private final Parameter precisionParameter;
        private final Parameter lambdaParameter;
        private final Parameter smoothedParameter;
        private final Parameter betaParameter;
        private final MatrixParameter dMatrix;
        protected SymmTridiagMatrix weightMatrix;
        protected SymmTridiagMatrix storedWeightMatrix;
    }

    class Covariates extends Default {
        public static final boolean NEW_APPROACH = true;

        Covariates( Parameter smoothingStatistic, Parameter precisionParameter, Parameter lambdaParameter,
                   List<Parameter> beta, MatrixParameter dMatrix, List<MatrixParameter> covariates, Parameter betaParameter, List<Parameter> deltaList) {
            super( smoothingStatistic, precisionParameter, lambdaParameter, betaParameter, dMatrix);

            if (beta != null) {
                for (Parameter betaParam : beta) {
                    addVariable(betaParam);
                }
            }

            if (deltaList != null) {
                this.delta = deltaList;
            } else {
                delta = new ArrayList<Parameter>();
                if (beta != null) {
                    for (int i = 0; i < beta.size(); i++) {
                        Parameter deltaParam = new Parameter.Default(1.0);
                        deltaParam.setParameterValue(0, 1.0);
                        delta.add(deltaParam);
                    }
                }
            }

            if (covariates != null) {
                for (MatrixParameter cov : covariates) {
                    addVariable(cov);
                }
            }
            this.beta = beta;
            this.covariates = covariates;
        }

        public List<Parameter> getBetaListParameter() {
            return this.beta;
        }

        public List<MatrixParameter> getCovariates() {
            return this.covariates;
        }

        protected void updateGammaWithCovariates(DenseVector currentGamma) {

            assert (beta != null);

            // Handle betaParameter / designMatrix

            if (NEW_APPROACH) {

                final int N = currentGamma.size();
                double[] update = new double[N];

                if (dMatrix != null) {
                    final int K = dMatrix.getColumnDimension();

                    if (N != dMatrix.getRowDimension()) {
                        throw new RuntimeException("Incorrect covariate dimensions (" + N + " != "
                                + dMatrix.getRowDimension() + ")");
                    }

                    for (int i = 0; i < N; ++i) {
                        for (int j = 0; j < K; ++j) {
                            update[i] += dMatrix.getParameterValue(i, j) * betaParameter.getParameterValue(j);
                        }
                    }
                }

                if (covariates != null) {
                    if (beta.size() != covariates.size()) {
                        throw new RuntimeException("beta.size(" + beta.size() + ") != covariates.size(" + covariates.size() + ")");
                    }

                    for (int k = 0; k < beta.size(); ++k) {

                        Parameter b = beta.get(k);
                        Parameter d = delta.get(k);
                        final int J = b.getDimension();
                        MatrixParameter covariate = covariates.get(k);
                        boolean transposed = isTransposed(N, J, covariate);

                        for (int i = 0; i < N; ++i) {
                            for (int j = 0; j < J; ++j) {
                                update[i] += covariate.getParameterValue(j, i) * b.getParameterValue(j) * d.getParameterValue(j);
                            }
                        }
                    }
                }

                for (int i = 0; i < N; ++i) {
                    currentGamma.set(i, currentGamma.get(i) - update[i]);
                }

            } else {
                DenseVector currentBeta = new DenseVector(beta.size());

                for (int i = 0; i < beta.size(); i++) {
                    currentBeta.set(i, beta.get(i).getParameterValue(0) * delta.get(i).getParameterValue(0));
                }

                //int numMissing = fieldLength - lastObservedIndex;
                //DenseVector tempVectCov = new DenseVector(numMissing);

                //System.err.println("covariates.size(): " + covariates.size());
                //System.err.println("covariates.get(0).getColumnDimension: " + covariates.get(0).getColumnDimension());
                //System.err.println("covariates.get(0).getRowDimension: " + covariates.get(0).getRowDimension());

                if (covariates != null) {

                    for (int i = 0; i < covariates.size(); i++) {
                        for (int j = 0; j < covariates.get(i).getColumnDimension(); j++) {
                            // System.err.println("j: " + j);
                            // System.err.println("covariates.get(i).getParameterValue(0,j): " + covariates.get(i).getParameterValue(0,j));
                            currentGamma.set(j, currentGamma.get(j) - covariates.get(i).getParameterValue(0, j) * currentBeta.get(i));
                        }
                    }
                }

                throw new RuntimeException("Should not get here.");
            }
        }

        public double[] getGradientWrtRegressionCoefficients() {


            if (beta == null) return null;

            if (beta.size() > 1 || getCovariates().size() > 1) {
                throw new RuntimeException("This is not the way to handle multidimensional parameters");
            }
            int numGridPoints = fieldLength - 1;

            Parameter b = beta.get(0);
            MatrixParameter covk = getCovariates().get(0);
            boolean transposed = isTransposed(numGridPoints + 1, b.getDimension(), covk);

            // TODO I believe we need one Parameter beta (is the same across loci) and List covariates (can differ across loci)

            double[] gradLogDens = new double[b.getDimension()];
            double[] gamma = getMeanAdjustedGamma().getData();

            double currentPrec = getPrecisionParameter().getParameterValue(0);

            for (int k = 0; k < b.getDimension(); k++) {

                double gradient = // numGridPoints / 2
                        +currentPrec * (gamma[0] - gamma[1]) * getCovariateValue(covk, 0, k, transposed) //covk.getParameterValue(k,0)
                                + currentPrec * (gamma[numGridPoints] - gamma[numGridPoints - 1]) * getCovariateValue(covk, numGridPoints, k, transposed) //covk.getParameterValue(k, numGridPoints)
//                    - 0.5 * currentPrec * (gamma[1] - gamma[0]) * (gamma[1] - gamma[0])
                        ;

                for (int i = 1; i < numGridPoints; i++) {
                    gradient +=
//                        -0.5 * currentPrec * (gamma[i + 1] - gamma[i]) * (gamma[i + 1] - gamma[i])
                            +currentPrec * (-gamma[i - 1] + 2 * gamma[i] - gamma[i + 1]) * getCovariateValue(covk, i, k, transposed); //covk.getParameterValue(k, i);
                }

                gradLogDens[k] = gradient;
            }

            return gradLogDens;
        }

        public double[] getDiagonalHessianWrtRegressionCoefficients() {

            if (beta == null) return null;
            int numGridPoints = fieldLength - 1;
            if (beta.size() > 1 || getCovariates().size() > 1) {
                throw new RuntimeException("This is not the way to handle multidimensional parameters");
            }

            Parameter b = beta.get(0);
            MatrixParameter covk = getCovariates().get(0);
            boolean transposed = isTransposed(numGridPoints + 1, b.getDimension(), covk);

            double[] hessian = new double[b.getDimension()];
            double currentPrec = getPrecisionParameter().getParameterValue(0);

            for (int k = 0; k < b.getDimension(); k++) {

                double h = // numGridPoints / 2
                        -currentPrec * (
                                getCovariateValue(covk, 0, k, transposed) //covk.getParameterValue(k, 0)
                                        - getCovariateValue(covk, 1, k, transposed)) //covk.getParameterValue(k,1))
                                * getCovariateValue(covk, 0, k, transposed) //covk.getParameterValue(k,0)
                                - currentPrec * (
                                getCovariateValue(covk, numGridPoints, k, transposed) //covk.getParameterValue(k, numGridPoints)
                                        - getCovariateValue(covk, numGridPoints - 1, k, transposed)) //covk.getParameterValue(k, numGridPoints - 1))
                                * getCovariateValue(covk, numGridPoints, k, transposed) //covk.getParameterValue(k, numGridPoints)
                        ;

                for (int i = 1; i < numGridPoints; i++) {
                    h += -currentPrec * (
                            -getCovariateValue(covk, i - 1, k, transposed) //covk.getParameterValue(k, i - 1)
                                    + 2 * getCovariateValue(covk, i, k, transposed) //covk.getParameterValue(k, i)
                                    - getCovariateValue(covk, i + 1, k, transposed)) //covk.getParameterValue(k, i + 1))
                            * getCovariateValue(covk, i, k, transposed) //covk.getParameterValue(k, i)
                    ;

                }

                hessian[k] = h;
            }

            return hessian;
        }


        private double getCovariateValue(MatrixParameter matrix, int i, int j, boolean transposed) {
            if (transposed) {
                return matrix.getParameterValue(j, i);
            } else {
                return matrix.getParameterValue(i, j);
            }
        }

        private boolean isTransposed(int N, int J, MatrixParameter matrix) {
            if (J == matrix.getRowDimension() && N == matrix.getColumnDimension()) {
                return true;
            } else if (J == matrix.getColumnDimension() && N == matrix.getRowDimension()) {
                return false;
            } else {
                throw new RuntimeException("Incorrect dimensions in " + matrix.getId() + " (r=" + matrix.getRowDimension() +
                        ",c=" + matrix.getColumnDimension() + ")");
            }
        }

        protected MatrixParameter dMatrix;
        final List<MatrixParameter> covariates;
        private final List<Parameter> beta;
        protected Parameter betaParameter;
        private final List<Parameter> delta;

    }

    class MissingCovariates extends Covariates {

        MissingCovariates(
                          Parameter smoothingStatistic,
                          Parameter precisionParameter,
                          Parameter lambdaParameter,
                          List<Parameter> beta,
                          MatrixParameter dMatrix,
                          List<MatrixParameter> covariates,
                          Parameter betaParameter,
                          List<Parameter> delta,
                          List<Parameter> covPrecParametersRecent,
                          List<Parameter> covPrecParametersDistant,
                          List<Parameter> firstObservedIndexParameter,
                          List<Parameter> lastObservedIndexParameter,
                          Parameter recentIndices,
                          Parameter distantIndices) {
            super( smoothingStatistic, precisionParameter, lambdaParameter, beta, dMatrix, covariates, betaParameter, delta);

            this.covPrecParametersRecent = covPrecParametersRecent;
            if (covPrecParametersRecent != null) {
                for (Parameter covPrecRecent : covPrecParametersRecent) {
                    addVariable(covPrecRecent);
                }
            }
            this.covPrecParametersDistant = covPrecParametersDistant;
            if (covPrecParametersDistant != null) {
                for (Parameter covPrecDistant : covPrecParametersDistant) {
                    addVariable(covPrecDistant);
                }
            }

            if (firstObservedIndexParameter != null) {
                this.firstObservedIndex = new int[firstObservedIndexParameter.size()];
                for (int i = 0; i < firstObservedIndexParameter.size(); i++) {
                    this.firstObservedIndex[i] = (int) firstObservedIndexParameter.get(i).getParameterValue(0);
                }

                if (recentIndices != null) {
                    // indices specify which covariates require default unobserved covariate data prior
                    this.recIndices = new int[firstObservedIndexParameter.size()];
                    for (int i = 0; i < firstObservedIndexParameter.size(); i++) {
                        this.recIndices[i] = (int) recentIndices.getParameterValue(i);
                    }
                } else {
                    // If specific covariates not specified by indices, need default unobserved covariate data prior for all covariates
                    this.recIndices = new int[firstObservedIndexParameter.size()];
                    for (int i = 0; i < firstObservedIndexParameter.size(); i++) {
                        this.recIndices[i] = i + 1;
                    }
                }

            }

            if (lastObservedIndexParameter != null) {
                this.lastObservedIndex = new int[lastObservedIndexParameter.size()];
                for (int i = 0; i < lastObservedIndexParameter.size(); i++) {
                    this.lastObservedIndex[i] = (int) lastObservedIndexParameter.get(i).getParameterValue(0);
                }

                if (distantIndices != null) {
                    // indices specify which covariates require default unobserved covariate data prior
                    this.distIndices = new int[lastObservedIndexParameter.size()];
                    for (int i = 0; i < lastObservedIndexParameter.size(); i++) {
                        this.distIndices[i] = (int) distantIndices.getParameterValue(i);
                    }
                } else {
                    // If specific covariates not specified by indices, need default unobserved covariate data prior for all covariates
                    this.distIndices = new int[lastObservedIndexParameter.size()];
                    for (int i = 0; i < lastObservedIndexParameter.size(); i++) {
                        this.distIndices[i] = i + 1;
                    }
                }

            }
            setupGMRFWeightsForMissingCov();
        }

        @Override
        protected double handleMissingValues() {

            assert (covPrecParametersRecent != null);
            assert (covariates != null);
            assert (covPrecParametersDistant != null);

            int numMissing;
            DenseVector tempVectMissingCov;
            SymmTridiagMatrix missingCovQ;
            DenseVector tempVectMissingCov2;
            int numMissingRecent;

            double currentLike = 0.0;

            if (lastObservedIndex != null) {
                for (int i = 0; i < covPrecParametersDistant.size(); i++) {

                    numMissing = fieldLength - lastObservedIndex[i];
                    tempVectMissingCov = new DenseVector(numMissing);
                    tempVectMissingCov2 = new DenseVector(numMissing);

                    missingCovQ = getScaledWeightMatrixForMissingCovDistant(covPrecParametersDistant.get(i).getParameterValue(0), i,
                            lastObservedIndex[i]);

                    for (int j = 0; j < numMissing; j++) {
                        tempVectMissingCov.set(j, covariates.get(distIndices[i] - 1).getParameterValue(0, lastObservedIndex[i] + j) -
                                covariates.get(distIndices[i] - 1).getParameterValue(0, lastObservedIndex[i] - 1));
                    }

                    missingCovQ.mult(tempVectMissingCov, tempVectMissingCov2);
                    currentLike += 0.5 * (numMissing) * Math.log(covPrecParametersDistant.get(i).getParameterValue(0))
                            - 0.5 * tempVectMissingCov.dot(tempVectMissingCov2);
                }
            }

            if (firstObservedIndex != null) {

                for (int i = 0; i < covPrecParametersRecent.size(); i++) {

                    numMissingRecent = firstObservedIndex[i] - 1;
                    tempVectMissingCov = new DenseVector(numMissingRecent);
                    tempVectMissingCov2 = new DenseVector(numMissingRecent);

                    missingCovQ = getScaledWeightMatrixForMissingCovRecent(covPrecParametersRecent.get(i).getParameterValue(0), i,
                            firstObservedIndex[i]);

                    for (int j = 0; j < numMissingRecent; j++) {
                        tempVectMissingCov.set(j, covariates.get(recIndices[i] - 1).getParameterValue(0, j) -
                                covariates.get(recIndices[i] - 1).getParameterValue(0, firstObservedIndex[i] - 1));
                    }

                    missingCovQ.mult(tempVectMissingCov, tempVectMissingCov2);
                    currentLike += 0.5 * (numMissingRecent) * Math.log(covPrecParametersRecent.get(i).getParameterValue(0))
                            - 0.5 * tempVectMissingCov.dot(tempVectMissingCov2);
                }
            }
            return currentLike;
        }

        private void setupGMRFWeightsForMissingCov() {

            if (firstObservedIndex != null) {
                weightMatricesForMissingCovRecent = new ArrayList<>();

                for (int i = 0; i < covPrecParametersRecent.size(); i++) {
                    double[] offdiagRec = new double[firstObservedIndex[i] - 2];
                    double[] diagRec = new double[firstObservedIndex[i] - 1];

                    for (int k = 0; k < firstObservedIndex[i] - 2; k++) {
                        offdiagRec[k] = -1;
                    }

                    for (int k = 1; k < firstObservedIndex[i] - 1; k++) {
                        diagRec[k] = 2.0;
                    }
                    diagRec[0] = 1.0;

                    weightMatricesForMissingCovRecent.add(i, new SymmTridiagMatrix(diagRec, offdiagRec));
                }

            }

            if (lastObservedIndex != null) {
                weightMatricesForMissingCovDistant = new ArrayList<>();

                for (int i = 0; i < covPrecParametersDistant.size(); i++) {
                    double[] offdiag = new double[fieldLength - lastObservedIndex[i] - 1];
                    double[] diag = new double[fieldLength - lastObservedIndex[i]];

                    //First set up the offdiagonal entries;

                    for (int k = 0; k < fieldLength - lastObservedIndex[i] - 1; k++) {
                        offdiag[k] = -1;
                    }

                    //Then set up the diagonal entries;
                    for (int k = 0; k < fieldLength - lastObservedIndex[i] - 1; k++) {
                        //	diag[i] = -(offdiag[i] + offdiag[i - 1]);
                        diag[k] = 2.0;
                    }
                    //Take care of the endpoint
                    diag[fieldLength - lastObservedIndex[i] - 1] = 1.0;

                    weightMatricesForMissingCovDistant.add(i, new SymmTridiagMatrix(diag, offdiag));
                }
            }

        }

        private SymmTridiagMatrix getScaledWeightMatrixForMissingCovRecent(double precision, int covIndex, int firstObs) {
            SymmTridiagMatrix a = weightMatricesForMissingCovRecent.get(covIndex).copy();
            for (int i = 0; i < a.numRows() - 1; i++) {
                a.set(i, i, a.get(i, i) * precision);
                a.set(i + 1, i, a.get(i + 1, i) * precision);
            }
            a.set(firstObs - 2, firstObs - 2,
                    a.get(firstObs - 2, firstObs - 2) * precision);
            return a;
        }

        private SymmTridiagMatrix getScaledWeightMatrixForMissingCovDistant(double precision, int covIndex, int lastObs) {
            SymmTridiagMatrix a = weightMatricesForMissingCovDistant.get(covIndex).copy();
            for (int i = 0; i < a.numRows() - 1; i++) {
                a.set(i, i, a.get(i, i) * precision);
                a.set(i + 1, i, a.get(i + 1, i) * precision);
            }
            a.set(fieldLength - lastObs - 1, fieldLength - lastObs - 1,
                    a.get(fieldLength - lastObs - 1, fieldLength - lastObs - 1) * precision);
            return a;
        }

        private final List<Parameter> covPrecParametersRecent;
        private final List<Parameter> covPrecParametersDistant;

        private List<SymmTridiagMatrix> weightMatricesForMissingCovRecent;
        private List<SymmTridiagMatrix> weightMatricesForMissingCovDistant;

        private int[] firstObservedIndex;
        private int[] lastObservedIndex;
        private int[] recIndices;
        private int[] distIndices;

    }

}
