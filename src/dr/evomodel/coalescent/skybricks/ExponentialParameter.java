package dr.evomodel.coalescent.skybricks;


import dr.inference.model.*;

public class ExponentialParameter extends Parameter.Abstract implements VariableListener {
    private final Parameter p;
    private final Bounds<Double> bounds;
    private int dimension;

    public ExponentialParameter(Parameter p) {
        this.p = p;
        p.addVariableListener(this);
        dimension = p.getDimension();

        bounds  = new DefaultBounds(Double.POSITIVE_INFINITY,2.0,p.getDimension());
    }

    @Override
    public boolean isImmutable() {
        return true; //I think
    }

    public int getDimension() {
        return dimension;
    }

    protected void storeValues() {
        p.storeParameterValues();

    }

    protected void restoreValues() {
        p.restoreVariableValues();
    }

    protected void acceptValues() {
        p.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }

    public double getParameterValue(int dim) {
        return Math.exp(p.getParameterValue(dim));

    }

    public void setParameterValue(int dim, double value) {
        throw new UnsupportedOperationException("Not implemented");
    }

    public void setParameterValueQuietly(int dim, double value) {
        throw new UnsupportedOperationException("Not implemented");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        throw new UnsupportedOperationException("Not implemented");
    }

    public String getParameterName() {
        if (getId() == null) {
            StringBuilder sb = new StringBuilder("exp");
            sb.append(".").append(p.getId());
            setId(sb.toString());
        }
        return getId();
    }

    public void addBounds(Bounds bounds) {
        throw new UnsupportedOperationException("Not implemented");
    }

    public Bounds<Double> getBounds() {
        return bounds;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent();
    }

}
