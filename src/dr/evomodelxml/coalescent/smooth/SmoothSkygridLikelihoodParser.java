package dr.evomodelxml.coalescent.smooth;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.smooth.SmoothSkygridLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class SmoothSkygridLikelihoodParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "smoothSkygridLikelihood";
    private static final String POPULATION_PARAMETER = GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER;
    private static final String PRECISION_PARAMETER = GMRFSkyrideLikelihoodParser.PRECISION_PARAMETER;
    private static final String INTERVALS = GMRFSkyrideLikelihoodParser.INTERVALS;
    private static final String POPULATION_TREE = GMRFSkyrideLikelihoodParser.POPULATION_TREE;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<IntervalList> intervalList = new ArrayList<>();
        if (xo.hasChildNamed(INTERVALS)) {
            XMLObject cxo = xo.getChild(INTERVALS);
            for (int i = 0; i < cxo.getChildCount(); ++i) {
                intervalList.add((IntervalList) cxo.getChild(i));
            }
        } else {
            XMLObject cxo = xo.getChild(INTERVALS);
            for (int i = 0; i < cxo.getChildCount(); ++i) {
                TreeModel tree = (TreeModel) xo.getChild(i);
                intervalList.add(new BigFastTreeIntervals(tree));
            }
        }

        return new SmoothSkygridLikelihood(null, null,
                null , null);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return SmoothSkygridLikelihood.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new XORRule(
                    new ElementRule(INTERVALS, new XMLSyntaxRule[]{
                            new ElementRule(IntervalList.class, 1, Integer.MAX_VALUE)
                    }),
                    new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                            new ElementRule(TreeModel.class, 1, Integer.MAX_VALUE)
                    })
            ),
    };
}
