package dr.evomodelxml.coalescent;


import dr.evomodel.coalescent.GMRFSkygrowthLikelihood;
import dr.xml.*;

public class GMRFSkygrowthGrowthRateStatisticParser  extends AbstractXMLObjectParser {

    public static final String SKYGROWTH_GROWTHRATE_STATISTIC = "GMRFSkygrowthGrowthRateStatistic";



    public String getParserName() { return SKYGROWTH_GROWTHRATE_STATISTIC; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GMRFSkygrowthLikelihood likelihood= (GMRFSkygrowthLikelihood) xo.getChild(GMRFSkygrowthLikelihood.class);


        return  likelihood.getGrowthRateStatistic();
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic for logging population sizes at set points in a skyline  ";
    }

    public Class getReturnType() { return GMRFSkygrowthGrowthRateStatisticParser.class; }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(GMRFSkygrowthLikelihood.class),
        };
    }
}
