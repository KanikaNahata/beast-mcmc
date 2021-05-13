package dr.evomodel.coalescent.skybricks;

import dr.inference.model.*;

public class SkyRateParameter extends Parameter.Abstract implements VariableListener, ModelListener {
    private final Parameter popSizes;
    private final boolean lastFlat;
    private final EpochProvider epochProvider;
    private final int dimension;
    private final Bounds<Double> bounds;

    public SkyRateParameter(Parameter popSizes, EpochProvider epochs) {
        this.popSizes = popSizes;
        popSizes.addVariableListener(this);

        epochProvider = epochs;
        // TODO figure out how to wrap all this up
        if (epochProvider instanceof Model) {
            ((Model) epochProvider).addModelListener(this);
        }
        lastFlat = epochProvider instanceof FixedGridEpochProvider;
        dimension = lastFlat? epochs.getEpochCount(): epochs.getEpochCount()-1;

        bounds  = new DefaultBounds(Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY,dimension);
    }

    @Override
    public boolean isImmutable() {
        return true; //I think
    }

    public int getDimension() {
        return dimension;
    }

    protected void storeValues() {
        popSizes.storeParameterValues();

    }

    protected void restoreValues() {
        popSizes.restoreVariableValues();
    }

    protected void acceptValues() {
        popSizes.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }

    public double getParameterValue(int dim) {
        if (dim == dimension-1 && lastFlat) {
            return 0;
        }
        return (popSizes.getParameterValue(dim) - popSizes.getParameterValue(dim + 1)) / epochProvider.getEpochDuration(dim);

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
            StringBuilder sb = new StringBuilder("rate");
                sb.append(".").append(popSizes.getId());
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

    /**
     * The model has changed. The model firing the event can optionally
     * supply a reference to an object and an index if appropriate. Use
     * of this extra information will be contingent on recognising what
     * model it was that fired the event.
     *
     * @param model
     * @param object
     * @param index
     */
    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        fireParameterChangedEvent();
    }

    /**
     * The model has been restored.
     * Required only for notification of non-models (say pure likelihoods) which depend on
     * models.
     *
     * @param model
     */
    @Override
    public void modelRestored(Model model) {
        //Do nothing I think
    }
}
