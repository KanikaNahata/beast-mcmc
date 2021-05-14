package dr.evomodelxml.coalescent.skybrick;

import dr.evomodel.coalescent.skybricks.ExponentialParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

public class ExponentialParameterParser extends AbstractXMLObjectParser {
    public static final String EXPONENTIAL_PARAMETER = "exponentialParameter";
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        return new ExponentialParameter(parameter);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
    };

    public String getParserDescription() {
        return "A parameter that exponentiates a wrapped parameter";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return EXPONENTIAL_PARAMETER;
    }
}
