package dr.evomodel.coalescent.skybricks;


import dr.inference.model.Parameter;

import java.util.Arrays;
import java.util.List;

public class FlexibleRootGridEpoch extends EpochProvider.Abstract {
    private double[] storedGridPoints;
    private final Parameter rootHeightParameter;

    private boolean storedEpochsKnown;
    public static final String FLEXIBLE_GRID_EPOCH = "flexibleGridEpoch";

    public FlexibleRootGridEpoch(Parameter rootHeightParameter,int numGridPoints, double cutOff){

        super(FLEXIBLE_GRID_EPOCH);
        epochCount = numGridPoints+1;

        gridPoints = new double[numGridPoints+1];
        storedGridPoints = new double[numGridPoints+1];

        Arrays.fill(gridPoints, 0);
        for (int pt = 0; pt < numGridPoints; pt++) {
            gridPoints[pt] = (pt) * (cutOff / numGridPoints);
        }
        gridPoints[numGridPoints] =rootHeightParameter.getParameterValue(0);
        epochsKnown=true;

        this.rootHeightParameter = rootHeightParameter;
        addVariable(rootHeightParameter);
    }

    public FlexibleRootGridEpoch(Parameter rootHeightParameter,Double[] grid){
        super("FLEXIBLE_GRID_EPOCH");

        List<Double> gridList =  Arrays.asList(grid);
        gridList.add(rootHeightParameter.getParameterValue(0));

        gridPoints =  new double[gridList.size()];
        storedGridPoints = new double[gridList.size()];
        for(int i = 0; i < gridList.size(); i++){
            gridPoints[i] = gridList.get(i);
        }
        epochCount = gridPoints.length;
        this.epochsKnown=true;
        this.rootHeightParameter=rootHeightParameter;
        addVariable(rootHeightParameter);
    }


    protected void calculateEpochs(){
        gridPoints[epochCount-1] = rootHeightParameter.getParameterValue(0);
        epochsKnown=true;
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {
        storedEpochsKnown=epochsKnown;
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

    protected void acceptState(){};

}
