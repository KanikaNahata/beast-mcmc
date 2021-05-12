package dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.xml.*;

public class RootHeightProxyParameter extends Parameter.Proxy {


    private final TreeModel tree;


    public RootHeightProxyParameter(String name,
                                    TreeModel tree) {
        super(name,  1);
        this.tree = tree;
    }


    @Override
    public double getParameterValue(int dim) {
        return tree.getNodeHeight(tree.getRoot());
    }

    @Override
    public void setParameterValue(int dim, double value) {
        if(dim>0){
            throw new IllegalArgumentException("Root parameter has dimension 1");
        }
        tree.setNodeHeight(tree.getRoot(), value);
        tree.pushTreeChangedEvent(tree.getRoot());
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        if(dim>0){
            throw new IllegalArgumentException("Root parameter has dimension 1");
        }
        tree.setNodeHeightQuietly(tree.getRoot(), value);
    }

    public String toString() {
        StringBuilder buffer = new StringBuilder(String.valueOf(getParameterValue(0)));
        Bounds bounds = null;

        for (int i = 1; i < getDimension(); i++) {
            buffer.append("\t").append(String.valueOf(getParameterValue(i)));
        }
        return buffer.toString();
    }

    @Override
    public void fireParameterChangedEvent() {
        tree.pushTreeChangedEvent();
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        setParameterValue(dim, value);
    }

    private static final String ROOT_HEIGHT_PARAMETER = "rootHeightProxyParameter";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
            return new dr.evomodel.treedatalikelihood.discrete.RootHeightProxyParameter(ROOT_HEIGHT_PARAMETER, tree);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class),
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter.class;
        }

        @Override
        public String getParserName() {
            return ROOT_HEIGHT_PARAMETER;
        }
    };
}