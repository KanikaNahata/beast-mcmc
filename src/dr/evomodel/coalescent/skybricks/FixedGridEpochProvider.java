package dr.evomodel.coalescent.skybricks;

import dr.inference.model.Parameter;

import java.util.Arrays;
import java.util.List;

public class FixedGridEpochProvider implements EpochProvider {

    private final double[] gridPoints;
    private final int epochCount;


    public FixedGridEpochProvider(int numGridPoints, double cutOff) {
        epochCount = numGridPoints+1;

        gridPoints = new double[numGridPoints+2];

        Arrays.fill(gridPoints, 0);
        for (int pt = 1; pt < numGridPoints+1; pt++) {
        gridPoints[pt] = (pt) * (cutOff / numGridPoints);
        }
        gridPoints[epochCount] =Double.POSITIVE_INFINITY;
    }

    public FixedGridEpochProvider(Double[] grid){
        boolean need0 = grid[0]!=0;
        boolean needInf = grid[grid.length-1]!=Double.POSITIVE_INFINITY;

        List<Double> gridList =  Arrays.asList(grid);
        if(need0){
            gridList.add(0,0.0);
        }
        if(needInf){
            gridList.add(Double.POSITIVE_INFINITY);
        }

        gridPoints =  new double[gridList.size()];
        for(int i = 0; i < gridList.size(); i++){
            gridPoints[i] = gridList.get(i);
        }
        epochCount = gridPoints.length-1;

    }

    /**
     * Get the epoch index corresponding to the time
     * @param time
     * @return
     */
    @Override
    public int getEpoch(double time){

            //https://stackoverflow.com/questions/6553970/find-the-first-element-in-a-sorted-array-that-is-greater-than-the-target
            //https://www.geeksforgeeks.org/first-strictly-greater-element-in-a-sorted-array-in-java/
            int low = 0, high = gridPoints.length-1;
            int ans = -1;
            while (low <= high) {

                int mid = (low + high) / 2;
                if (gridPoints[mid] <= time) {
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
                ans = epochCount;
            }
            return ans-1; //refers to the epoch that ends at this index

        }
    @Override
    public int getEpochCount(){
        return epochCount;
    };

    @Override
    public double getEpochStartTime(int i) {
        return gridPoints[i];
    }
    @Override
    public double getEpochEndTime(int i){
        return gridPoints[i + 1];
    }
    @Override
    public double getEpochDuration(int i){
        return gridPoints[i+1]-gridPoints[i];
    }


}
