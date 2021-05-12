package dr.evomodelxml.coalescent.skybrick;

import dr.evolution.coalescent.IntervalList;
import dr.evomodel.coalescent.skybricks.EpochProvider;
import dr.evomodel.coalescent.skybricks.SkybrickLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

public class SkybrickLikelihoodParser extends AbstractXMLObjectParser {
    public static final String SKYBRICK = "skybrickLikelihood";
    public static final String INTERVALS = "intervals";


    public static final String POPULATION_PARAMETER = "populationSizes";
    public static final String EPOCHS = "epochs";
    public static final String CONSTANT = "constant";
    public static final String EXPONENTIAL = "exponential";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject cxo = xo.getChild(POPULATION_PARAMETER);
        Parameter popParameter = (Parameter) cxo.getChild(Parameter.class);
        cxo = xo.getChild(EPOCHS);
        EpochProvider epochProvider = (EpochProvider) cxo.getChild(EpochProvider.class);

        IntervalList intervalList = (IntervalList) xo.getChild(INTERVALS).getChild(IntervalList.class);


        return new SkybrickLikelihood(epochProvider, popParameter, intervalList, xo.getAttribute("piecewiseApproximation", CONSTANT));
    }

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(EPOCHS, new XMLSyntaxRule[]{
                new ElementRule(EpochProvider.class)
                }),

                new ElementRule(INTERVALS, new XMLSyntaxRule[]{
                        new ElementRule(IntervalList.class)
                }, "The interval list from which the coalescent likelihood will be calculated")
        };
    }

    @Override
    public String getParserDescription() {
        return "Parses the a flexible piecewise coalescent";
    }

    @Override
    public Class getReturnType() {
        return SkybrickLikelihood.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return SKYBRICK;
    }
}
