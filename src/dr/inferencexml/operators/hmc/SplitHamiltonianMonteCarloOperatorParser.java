package dr.inferencexml.operators.hmc;

import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.SplitHamiltonianMonteCarloOperator;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 */

public class SplitHamiltonianMonteCarloOperatorParser extends AbstractXMLObjectParser { //todo: merge with HMC parser?

    private final static String SPLIT_HMC = "splitHamiltonianMonteCarloOperator";
    private final static String N_STEPS = "nSteps";
    private final static String N_INNER_STEPS = "nInnerSteps";
    private final static String STEP_SIZE = "stepSize";
    private final static String RELATIVE_SCALE = "relativeScale";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double stepSize = xo.getDoubleAttribute(STEP_SIZE);
        double relativeScale = xo.getDoubleAttribute(RELATIVE_SCALE);

        ReversibleHMCProvider reversibleHMCproviderA = (ReversibleHMCProvider) xo.getChild(0);
        ReversibleHMCProvider reversibleHMCproviderB = (ReversibleHMCProvider) xo.getChild(1);

        int nStep = xo.getAttribute(N_STEPS, 5);
        int nInnerStep = xo.getAttribute(N_INNER_STEPS, 5);

        return new SplitHamiltonianMonteCarloOperator(weight, reversibleHMCproviderA, reversibleHMCproviderB,
                stepSize, relativeScale, nStep
                , nInnerStep);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    protected final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(N_STEPS, true),
            AttributeRule.newIntegerRule(N_INNER_STEPS, true),
            AttributeRule.newDoubleRule(STEP_SIZE),
            AttributeRule.newDoubleRule(RELATIVE_SCALE)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return null;
    }

    @Override
    public String getParserName() {
        return SPLIT_HMC;
    }
}