package dr.evomodelxml.coalescent;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.GMRFSkygrowthLikelihood;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class GMRFSkygrowthLikelihoodParser extends AbstractXMLObjectParser {
    public static final String SKYGROWTH = "skygrowthLikelihood";
    public static final String POPULATION_TREE = "populationTree";
    public static final String INTERVALS = "intervals";
    public static final String INCLUDE = "include";
    public static final String EXCLUDE = "exclude";

    public static final String POPULATION_PARAMETER = "populationSizes";
    public static final String PRECISION_PARAMETER = "precisionParameter";
    public static final String GRID_POINTS = "gridPoints";
    public static final String NUM_GRID_POINTS = "numGridPoints";
    public static final String CUT_OFF = "cutOff";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        IntervalList intervalList = null;
        TreeModel treeModel = null;
        XMLObject cxo = xo.getChild(POPULATION_PARAMETER);
        Parameter popParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PRECISION_PARAMETER);
        Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

        Parameter gridPoints = null;
        if (xo.getChild(GRID_POINTS) != null) {
            cxo = xo.getChild(GRID_POINTS);
            gridPoints = (Parameter) cxo.getChild(Parameter.class);
        }



        if (xo.getChild(NUM_GRID_POINTS) != null) {
            cxo = xo.getChild(NUM_GRID_POINTS);
            Parameter numGridPoints = (Parameter) cxo.getChild(Parameter.class);
            cxo = xo.getChild(CUT_OFF);
            Parameter cutOff = (Parameter) cxo.getChild(Parameter.class);
            double step = cutOff.getParameterValue(0) / (numGridPoints.getParameterValue(0));
            double currentPos = 0;
            gridPoints = new Parameter.Default((int)numGridPoints.getParameterValue(0));
            for (int i = 0; i < gridPoints.getDimension(); i++) {
                currentPos += step;
                gridPoints.setParameterValue(i, currentPos);
            }
        }
        cxo = xo.getChild(POPULATION_TREE);
        if (cxo != null) {
            treeModel = (TreeModel) cxo.getChild(TreeModel.class);

            TaxonList includeSubtree = null;

            if (xo.hasChildNamed(INCLUDE)) {
                includeSubtree = (TaxonList) xo.getElementFirstChild(INCLUDE);
            }

            List<TaxonList> excludeSubtrees = new ArrayList<TaxonList>();

            if (xo.hasChildNamed(EXCLUDE)) {
                cxo = xo.getChild(EXCLUDE);
                for (int i = 0; i < cxo.getChildCount(); i++) {
                    excludeSubtrees.add((TaxonList) cxo.getChild(i));
                }
            }

            try {
                intervalList = new TreeIntervals(treeModel, includeSubtree, excludeSubtrees);
                // TreeIntervals now deals with all the interval stuff
//                return new CoalescentLikelihood(treeModel, includeSubtree, excludeSubtrees, demoModel);
            } catch (TreeUtils.MissingTaxonException mte) {
                throw new XMLParseException("treeModel missing a taxon from taxon list in " + getParserName() + " element");
            }


        } else {
            intervalList = (IntervalList) xo.getChild(INTERVALS).getChild(IntervalList.class);
        }

        return new GMRFSkygrowthLikelihood(gridPoints, precParameter, popParameter, intervalList);
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
                new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new XORRule(
                        new ElementRule(CoalescentLikelihoodParser.POPULATION_TREE, new XMLSyntaxRule[]{
                                new ElementRule(TreeModel.class)
                        }),
                        new ElementRule(INTERVALS, new XMLSyntaxRule[]{
                                new ElementRule(IntervalList.class)
                        }, "The interval list from which the coalescent likelihood will be calculated")

                ),
                new XORRule(
                        new ElementRule(GRID_POINTS, new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)
                        }),
                        new AndRule(
                                new ElementRule(NUM_GRID_POINTS, new XMLSyntaxRule[]{
                                        new ElementRule(Parameter.class)
                                }, "The interval list from which the coalescent likelihood will be calculated"),

                                new ElementRule(CUT_OFF, new XMLSyntaxRule[]{
                                        new ElementRule(Parameter.class),
                                }, "Cut off. ie. the last point in the past where the population is allowed to change"
                                )
                        )
                ),
        };
    }

    @Override
    public String getParserDescription() {
        return "Parses the a sky growth likelihood as described in Didelot and Volz 2017";
    }

    @Override
    public Class getReturnType() {
        return GMRFSkygrowthLikelihood.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    @Override
    public String getParserName() {
        return SKYGROWTH;
    }
}
