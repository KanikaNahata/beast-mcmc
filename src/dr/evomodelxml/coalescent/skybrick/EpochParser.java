package dr.evomodelxml.coalescent.skybrick;

import dr.evolution.coalescent.IntervalList;
import dr.evomodel.coalescent.skybricks.*;
import dr.inference.model.Parameter;
import dr.xml.*;


public class EpochParser extends AbstractXMLObjectParser {
    public static final String GROUP_SIZES = "groupSizes";
    public static final String INTERVALS = "intervals";
    public static final String GRID_POINTS = "gridPoints";
    public static final String FIXED_EPOCH = "fixedEpoch";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";
    public static final String EPOCH_MODEL = "epochModel";
    public static final String ROOT_HEIGHT_PARAMETER = "rootHeightParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // false no implemented yet but would be like skygrid but with root behind cutoff and last epoch defined as cutoff to root
//        boolean fixedEpoch = xo.getAttribute(FIXED_EPOCH, true);

        XMLObject cxo = xo.getChild(GROUP_SIZES);
        Parameter groupParameter = null;
        if (cxo != null) {
            groupParameter = (Parameter) cxo.getChild(Parameter.class);
        }
        IntervalList intervalList = null;
        if (xo.getChild(INTERVALS) != null) {
            cxo = xo.getChild(INTERVALS);
            intervalList = (IntervalList) cxo.getChild(IntervalList.class);
        }

        Parameter gridPoints = null;
        if (xo.getChild(GRID_POINTS) != null) {
            cxo = xo.getChild(GRID_POINTS);
            gridPoints = (Parameter) cxo.getChild(Parameter.class);
        }

        Parameter numGridPoints = null;
        if (xo.getChild(NUM_GRID_POINTS) != null) {
            cxo = xo.getChild(NUM_GRID_POINTS);
            numGridPoints = (Parameter) cxo.getChild(Parameter.class);
        }

        Parameter cutOff = null;
        if (xo.getChild(CUT_OFF) != null) {
            cxo = xo.getChild(CUT_OFF);
            cutOff = (Parameter) cxo.getChild(Parameter.class);
        }

        Parameter rootHeightParameter = null;
        if (xo.getChild(ROOT_HEIGHT_PARAMETER) != null) {
            cxo = xo.getChild(ROOT_HEIGHT_PARAMETER);
            rootHeightParameter = (Parameter) cxo.getChild(Parameter.class);
        }

        EpochProvider epochProvider = null;

        if (gridPoints != null) {
            if (rootHeightParameter != null) {
                epochProvider = new FlexibleRootGridEpoch(rootHeightParameter, gridPoints.getValues());
            } else {
                epochProvider = new FixedGridEpochProvider(gridPoints.getValues());
            }
        } else if (cutOff != null && numGridPoints != null) {
            if (rootHeightParameter != null) {
                epochProvider = new FlexibleRootGridEpoch(rootHeightParameter, numGridPoints.getValue(0).intValue(), cutOff.getValue(0));
            } else {
                epochProvider = new FixedGridEpochProvider(numGridPoints.getValue(0).intValue(), cutOff.getValue(0));
            }
        } else if (intervalList != null) {
            epochProvider = new SkylineEpochProvider(intervalList, groupParameter);
        }
        return epochProvider;
    }

    /**
     * @return an array of syntax rules required by this element.
     * Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new XORRule(
                        new XORRule(
                                new ElementRule(GRID_POINTS, new XMLSyntaxRule[]{
                                        new ElementRule(Parameter.class)
                                }),
                                new AndRule(new ElementRule(CUT_OFF, new XMLSyntaxRule[]{
                                        new ElementRule(Parameter.class)
                                }),
                                        new ElementRule(NUM_GRID_POINTS, new XMLSyntaxRule[]{
                                                new ElementRule(Parameter.class)
                                        }))
                        ),
                        new AndRule(new ElementRule(INTERVALS, new XMLSyntaxRule[]{
                                new ElementRule(IntervalList.class)
                        }),
                                new ElementRule(GROUP_SIZES, new XMLSyntaxRule[]{
                                        new ElementRule(Parameter.class)
                                })
                        )
                )
        };
    }

    @Override
    public String getParserDescription() {
        return "Parses the a flexible epoch model";
    }

    @Override
    public Class getReturnType() {
        return EpochProvider.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return EPOCH_MODEL;
    }
}

