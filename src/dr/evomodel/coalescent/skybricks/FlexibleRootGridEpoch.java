package dr.evomodel.coalescent.skybricks;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;
import java.util.List;

public class FlexibleRootGridEpoch extends AbstractModel implements EpochProvider {
    private double[] gridPoints;
    private double[] storedGridPoints;

    private final int epochCount;
    private final Parameter rootHeightParameter;
    private boolean rootHeightKnown;
    private boolean storedRootHieghtKnown;
    public FlexibleRootGridEpoch(Parameter rootHeightParameter,int numGridPoints, double cutOff){
        super("FlexibleRootGridEpoch");
        epochCount = numGridPoints+1;
        this.rootHeightParameter = rootHeightParameter;
        addVariable(rootHeightParameter);

        gridPoints = new double[numGridPoints+2];
        storedGridPoints = new double[numGridPoints+2];

        Arrays.fill(gridPoints, 0);
        for (int pt = 1; pt < numGridPoints+1; pt++) {
            gridPoints[pt] = (pt) * (cutOff / numGridPoints);
        }
        gridPoints[epochCount] =rootHeightParameter.getParameterValue(0);
        rootHeightKnown=true;
    }
    public FlexibleRootGridEpoch(Parameter rootHeightParameter,Double[] grid){
        super("FlexibleRootGridEpoch");
        boolean need0 = grid[0]!=0;


        List<Double> gridList =  Arrays.asList(grid);
        if(need0){
            gridList.add(0,0.0);
        }
        gridList.add(rootHeightParameter.getParameterValue(0));


        gridPoints =  new double[gridList.size()];
        storedGridPoints = new double[gridList.size()];
        for(int i = 0; i < gridList.size(); i++){
            gridPoints[i] = gridList.get(i);
        }
        epochCount = gridPoints.length-1;
        this.rootHeightKnown=true;
        this.rootHeightParameter=rootHeightParameter;
    }


    @Override
    public int getEpoch(double time){
        if(!rootHeightKnown){
            calculateGridPoints();
        }
        //https://stackoverflow.com/questions/6553970/find-the-first-element-in-a-sorted-array-that-is-greater-than-the-target
        //https://www.geeksforgeeks.org/first-strictly-greater-element-in-a-sorted-array-in-java/
        int low = 0, high = gridPoints.length-1;
        int ans = -1;
        while (low <= high) {

            int mid = (low + high) / 2;
            if (gridPoints[mid] < time) {
                /* This index, and everything below it, must not be the first element
                 * greater than what we're looking for because this element is no greater
                 * than the element.
                 */
                low = mid + 1;
            } else {

                ans = mid;
                high = mid - 1;
            }

        }
        if (ans == -1) {
            // There is no greater element so the "first greater"
            ans = epochCount;
        }
        if(ans==0){
            return ans;
        }
        return ans-1; //refers to the epoch that ends at this index

    }

    @Override
    public int getEpochCount(){
        return epochCount;
    };

    @Override
    public double getEpochStartTime(int i) {
        if(!rootHeightKnown){
            calculateGridPoints();
        }
        return gridPoints[i];
    }
    @Override
    public double getEpochEndTime(int i){
        if(!rootHeightKnown){
            calculateGridPoints();
        }
        return gridPoints[i + 1];
    }
    @Override
    public double getEpochDuration(int i){
        if(!rootHeightKnown){
            calculateGridPoints();
        }
        return gridPoints[i+1]-gridPoints[i];
    }

    private void calculateGridPoints(){
        gridPoints[epochCount] = rootHeightParameter.getParameterValue(0);
        rootHeightKnown=true;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {
        storedRootHieghtKnown=rootHeightKnown;
        System.arraycopy(gridPoints,0,storedGridPoints,0,gridPoints.length);
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {
        rootHeightKnown = storedRootHieghtKnown;
        double[] tmp1 = storedGridPoints;
        storedGridPoints = gridPoints;
        gridPoints = tmp1;
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {

    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     *
     * @param variable
     * @param index
     * @param type
     */
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        rootHeightKnown=false;
        fireModelChanged();
    }
}
