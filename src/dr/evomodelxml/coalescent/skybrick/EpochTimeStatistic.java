package dr.evomodelxml.coalescent.skybrick;

import dr.evomodel.coalescent.skybricks.EpochProvider;
import dr.inference.model.Parameter;
import dr.xml.*;

public class EpochTimeStatistic extends AbstractXMLObjectParser {
    public final static String EPOCH_TIME_STATISTIC = "epochTimeStatistic";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        EpochProvider epochs = (EpochProvider) xo.getChild(EpochProvider.class);
        return epochs.getStatistic("epochTime");
    }

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(EpochProvider.Abstract.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "returns a statistic that reports the transition times of an epoch model";
    }

    @Override
    public Class getReturnType() {
        return EpochProvider.Abstract.EpochTimes.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return EPOCH_TIME_STATISTIC;
    }
}
