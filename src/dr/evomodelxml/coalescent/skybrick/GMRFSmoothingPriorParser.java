package dr.evomodelxml.coalescent.skybrick;

import dr.evomodel.coalescent.skybricks.EpochProvider;
import dr.evomodel.coalescent.skybricks.GMRFSmoothingPrior;
import dr.inference.glm.GeneralizedLinearModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;
import java.util.List;

public class GMRFSmoothingPriorParser extends AbstractXMLObjectParser {
    public static final String GMRFSMOOTHING_PRIOR = "gmrfSmoothingPrior";
    public static final String SMOOTHED_PARAMETER = "smoothedParameter";
    public static final String PRECISION_PARAMETER = "precisionParameter";


    public static final String FIRST_OBSERVED_INDEX = "firstObservedIndex";
    public static final String LAST_OBSERVED_INDEX = "lastObservedIndex";
    public static final String COV_PREC_PARAM = "covariatePrecision";
    public static final String COV_PREC_REC = "covariatePrecisionRecent";
    public static final String COV_PREC_DIST = "covariatePrecisionDistant";
    public static final String REC_INDICES = "covIndicesMissingRecent";
    public static final String DIST_INDICES = "covIndicesMissingDistant";
    public static final String COVARIATES = "covariates";

    public static final String LAMBDA_PARAMETER = "lambdaParameter";
    public static final String BETA_PARAMETER = "betaParameter";
    public static final String DELTA_PARAMETER = "deltaParameter";
    public static final String SINGLE_BETA = "singleBeta";
    public static final String COVARIATE_MATRIX = "covariateMatrix";




    public static final String USE_GLM_MODEL = "useGlmModel";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        GMRFSmoothingPrior.Builder builder = new GMRFSmoothingPrior.Builder();

        XMLObject cxo = xo.getChild(SMOOTHED_PARAMETER);
        Parameter smoothedParameter = (Parameter) cxo.getChild(Parameter.class);
        builder.setSmoothingParameter(smoothedParameter);

        cxo = xo.getChild(PRECISION_PARAMETER);
        Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);
        builder.setPrecisionParameter(precParameter);

//TODO time aware smoothing is not currently implemented

//        EpochProvider epochProvider = (EpochProvider) xo.getChild(EpochProvider.class);
//        if (epochProvider != null) {
//            builder.setEpochModel(epochProvider);
//
//        }

        Parameter lambda;
        if (xo.getChild(LAMBDA_PARAMETER) != null) {
            cxo = xo.getChild(LAMBDA_PARAMETER);
            lambda = (Parameter) cxo.getChild(Parameter.class);
        } else {
            lambda = new Parameter.Default(LAMBDA_PARAMETER, 1.0);
        }
        builder.setLambdaParameter(lambda);
        List<Parameter> firstObservedIndex = null;
        if (xo.hasChildNamed(FIRST_OBSERVED_INDEX)) {
            firstObservedIndex  = new ArrayList<Parameter>();
            cxo = xo.getChild(FIRST_OBSERVED_INDEX);
            final int numInd = cxo.getChildCount();

            for(int i=0; i< numInd; ++i) {
                firstObservedIndex.add((Parameter) cxo.getChild(i));
            }
            builder.setFirstObservedIndexParameter(firstObservedIndex);
        }
        List<Parameter> lastObservedIndex = null;
        if (xo.hasChildNamed(LAST_OBSERVED_INDEX)) {
            lastObservedIndex  = new ArrayList<Parameter>();
            cxo = xo.getChild(LAST_OBSERVED_INDEX);
            final int numObsInd = cxo.getChildCount();

            for(int i=0; i< numObsInd; ++i) {
                lastObservedIndex.add((Parameter) cxo.getChild(i));
            }
            builder.setLastObservedIndexParameter(lastObservedIndex);
        }


        Parameter betaParameter = null;
        if (xo.hasChildNamed(SINGLE_BETA)) {
            betaParameter = (Parameter) xo.getElementFirstChild(SINGLE_BETA);
        }
        builder.setBetaParameter(betaParameter);

        List<Parameter> betaList = null;
        if (xo.getChild(BETA_PARAMETER) != null) {
            betaList = new ArrayList<Parameter>();
            cxo = xo.getChild(BETA_PARAMETER);
            final int numBeta = cxo.getChildCount();
            for (int i = 0; i < numBeta; ++i) {
                betaList.add((Parameter) cxo.getChild(i));
            }
        }
        builder.setBetaList(betaList);


        if (xo.getChild(DELTA_PARAMETER) != null) {
            List<Parameter> deltaList = new ArrayList<Parameter>();

            cxo = xo.getChild(DELTA_PARAMETER);
            final int numDelta = cxo.getChildCount();
            if (numDelta != betaList.size()) {
                throw new XMLParseException("Cannot have different number of delta and beta parameters");
            }
            for (int i = 0; i < numDelta; ++i) {
                deltaList.add((Parameter) cxo.getChild(i));
            }

            builder.setDeltaList(deltaList);
        }

        if (xo.getChild(COVARIATE_MATRIX) != null) {
            cxo = xo.getChild(COVARIATE_MATRIX);
            MatrixParameter   dMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);

            builder.setdMatrix(dMatrix);
        }
        List<Parameter>  covPrecParamRecent = new ArrayList<Parameter>();
        if(xo.hasChildNamed(COV_PREC_REC)){
            cxo = xo.getChild(COV_PREC_REC);
            for(int i = 0; i < cxo.getChildCount(); ++i){
                covPrecParamRecent.add((Parameter) cxo.getChild(i));
            }
            builder.setCovPrecParametersRecent(covPrecParamRecent);
        }


        List<Parameter>   covPrecParamDistant = new ArrayList<Parameter>();
        if(xo.hasChildNamed(COV_PREC_DIST)){
            cxo = xo.getChild(COV_PREC_DIST);
            for(int i = 0; i < cxo.getChildCount(); ++i){
                covPrecParamDistant.add((Parameter) cxo.getChild(i));
            }
            builder.setCovPrecParametersDistant(covPrecParamDistant);
        }

        if  (xo.hasChildNamed(COV_PREC_PARAM)){

            if(firstObservedIndex != null) {
                covPrecParamRecent = new ArrayList<Parameter>();
            }
            if(lastObservedIndex != null) {
                covPrecParamDistant = new ArrayList<Parameter>();
            }
            cxo = xo.getChild(COV_PREC_PARAM);

            for(int i=0; i < cxo.getChildCount(); ++i){
                if(firstObservedIndex != null) {
                    covPrecParamRecent.add((Parameter) cxo.getChild(i));
                }
                if(lastObservedIndex != null) {
                    covPrecParamDistant.add((Parameter) cxo.getChild(i));
                }
            }
            builder.setCovPrecParametersRecent(covPrecParamRecent);
            builder.setCovPrecParametersDistant(covPrecParamDistant);
        }

        if (xo.getChild(REC_INDICES) != null) {
            cxo = xo.getChild(REC_INDICES);
            Parameter  recentIndices = (Parameter) cxo.getChild(Parameter.class);
            builder.setRecentIndices(recentIndices);
        }

        if (xo.getChild(DIST_INDICES) != null) {
            cxo = xo.getChild(DIST_INDICES);
            Parameter  distantIndices = (Parameter) cxo.getChild(Parameter.class);
            builder.setDistantIndices(distantIndices);
        }

        if (xo.hasChildNamed(COVARIATES)){
            List<MatrixParameter>  covariates = new ArrayList<MatrixParameter>();
            cxo = xo.getChild(COVARIATES);
            final int numCov = cxo.getChildCount();

            for (int i = 0; i < numCov; ++i) {
                covariates.add((MatrixParameter) cxo.getChild(i));
            }
            builder.setCovariates(covariates);
        }


        GeneralizedLinearModel glm = (GeneralizedLinearModel) xo.getChild(GeneralizedLinearModel.class);
        builder.setGlm(glm);


        return builder.build();
    }

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return GMRFSmoothingPrior.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return GMRFSMOOTHING_PRIOR;
    }
}
