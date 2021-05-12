package dr.evomodelxml.coalescent.skybrick;

import dr.evomodel.coalescent.skybricks.EpochProvider;
import dr.evomodel.coalescent.skybricks.SkyRateParameter;
import dr.inference.model.Parameter;
import dr.xml.*;


public class SkyRateParameterParser extends AbstractXMLObjectParser {
    public static final String SKYRATE_PARAMETER = "skyrateParameter";
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        EpochProvider epochProvider = (EpochProvider) xo.getChild(EpochProvider.class);

        return new SkyRateParameter(parameter,epochProvider);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new ElementRule(EpochProvider.class)
    };

    public String getParserDescription() {
        return "A parameter that caluculates the rate of change between neighboring values in a wrapped parameter ";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return SKYRATE_PARAMETER;
    }
}


