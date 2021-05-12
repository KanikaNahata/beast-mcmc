package dr.evomodel.coalescent.skybricks;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evomodel.coalescent.OldAbstractCoalescentLikelihood;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;

public class SkylineEpochProvider extends AbstractModel implements EpochProvider {
    private final IntervalList intervalList;
    private final Parameter groupSizeParameter;
    private boolean epochsKnown;
    private boolean storedEpochsKnown;
    private double[] gridPoints;
    private double[] storedGridPoints;
    private final int epochCount;

    public SkylineEpochProvider(IntervalList intervalList, Parameter groupSizeParameter){
        super(GMRFSkyrideLikelihoodParser.SKYLINE_LIKELIHOOD); //For fun ;)
        this.intervalList = intervalList;
        if(intervalList instanceof Model){
            addModel((Model) intervalList);
        }

        this.groupSizeParameter = groupSizeParameter;
        addVariable(groupSizeParameter);

        epochCount = groupSizeParameter.getDimension();
        gridPoints = new double[groupSizeParameter.getDimension()+1];

    }

    @Override
    public int getEpoch(double time) {
        if(!epochsKnown){
            calculateGridPoints();
        }
        //https://stackoverflow.com/questions/6553970/find-the-first-element-in-a-sorted-array-that-is-greater-than-the-target
        //https://www.geeksforgeeks.org/first-strictly-greater-element-in-a-sorted-array-in-java/
        int low = 0, high = epochCount;
        int ans = -1;
        while (low <= high) {

            int mid = (low + high) / 2;
            if (getEpochEndTime(mid) <= time) {
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
            // There is no greater element so the "first greater" is outside the array.
            ans = epochCount + 1;
        }
        return ans;

    }

    @Override
    public int getEpochCount() {
        return epochCount;
    }

    @Override
    public double getEpochStartTime(int i) {
        if(!epochsKnown){
            calculateGridPoints();
        }
        return gridPoints[i];
    }

    @Override
    public double getEpochEndTime(int i) {
        if(!epochsKnown){
            calculateGridPoints();
        }
        return gridPoints[i+1];
    }

    @Override
    public double getEpochDuration(int i) {
        if(!epochsKnown){
            calculateGridPoints();
        }
        return gridPoints[i+1] - gridPoints[i] ;
    }

    private void calculateGridPoints(){
//         gridPoints = new double[groupSizeParameter.getDimension()+1];
        Arrays.fill(gridPoints,0.0);
        double timeEnd = 0.0;
        int groupIndex = 1;
        int subIndex = 0;
        for (int i = 0; i < intervalList.getIntervalCount(); i++) {

            timeEnd += intervalList.getInterval(i);

            if (intervalList.getIntervalType(i) == IntervalType.COALESCENT) {
                subIndex += 1;
                if (subIndex >= groupSizeParameter.getParameterValue(groupIndex)) {
                    gridPoints[groupIndex] = timeEnd;
                    groupIndex += 1;
                    subIndex = 0;
                }
            }
        }
        gridPoints[groupSizeParameter.getDimension()] = timeEnd;

    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        epochsKnown = false;
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {
        storedEpochsKnown = epochsKnown;
        System.arraycopy(gridPoints,0,storedGridPoints,0,gridPoints.length);
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {
        epochsKnown = storedEpochsKnown;
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
        epochsKnown=false;
    }
}
