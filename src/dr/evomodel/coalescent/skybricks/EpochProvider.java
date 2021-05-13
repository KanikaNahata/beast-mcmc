package dr.evomodel.coalescent.skybricks;

import dr.inference.model.*;

public interface EpochProvider extends StatisticList{
    int getEpoch(double time);

    int getEpochCount();

    double getEpochStartTime(int i);

    double getEpochEndTime(int i);

    double getEpochDuration(int i);



    abstract class   Abstract extends AbstractModel implements EpochProvider {
        protected double[] gridPoints = new double[0];
        protected int epochCount = 0;
        protected boolean epochsKnown;


        public Abstract(String name){
            super(name);
            addStatistic( new EpochTimes());
        }
        /**
         * A binary search epoch that ends in the first grid time equal to or greater than the specified
         * time
         * @param time
         * @return int
         */
        public final int getEpoch(double time){
            if(!epochsKnown){
                calculateEpochs();
            }
            //https://stackoverflow.com/questions/6553970/find-the-first-element-in-a-sorted-array-that-is-greater-than-the-target
            //https://www.geeksforgeeks.org/first-strictly-greater-element-in-a-sorted-array-in-java/
            // really getting the first gridpoint greater than or equal to a given time.

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
                }else if(gridPoints[mid]==time){
                    return mid;
                }else{
                    ans = mid;
                    high = mid - 1;
                }
            }
            return ans;
        }
        public int getEpochCount(){
            return epochCount;
        };

        public double getEpochStartTime(int i) {
            if(!epochsKnown){
                calculateEpochs();
            }
            if(i==0){
                return 0;
            }
            return gridPoints[i-1];
        }
        public double getEpochEndTime(int i){
            if(!epochsKnown){
                calculateEpochs();
            }
            return gridPoints[i];
        }
        public double getEpochDuration(int i){
            if(!epochsKnown){
                calculateEpochs();
            }
            if(i==0){
                return gridPoints[i];
            }
            return gridPoints[i]-gridPoints[i-1];
        }
        protected void calculateEpochs(){
            throw new UnsupportedOperationException("don't call the base class");
        };

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            epochsKnown =false;
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
            epochsKnown = false;
        }
        public class EpochTimes extends Statistic.Abstract {

            public EpochTimes() {
                super("epochTime");
            }

            public int getDimension() { return gridPoints.length; }

            public double getStatisticValue(int i) {
                return gridPoints[i];
            }
        }

    }
}
