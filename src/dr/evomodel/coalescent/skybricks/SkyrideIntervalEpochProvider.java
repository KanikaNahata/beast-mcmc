package dr.evomodel.coalescent.skybricks;

import dr.evolution.coalescent.IntervalList;

public class SkyrideIntervalEpochProvider implements EpochProvider {
    private final int epochCount;
    private final IntervalList intervals;

    public SkyrideIntervalEpochProvider(IntervalList intervals) {
        epochCount = intervals.getIntervalCount();
        this.intervals = intervals;
    }

    @Override
    public int getEpoch(double time) {
        //https://stackoverflow.com/questions/6553970/find-the-first-element-in-a-sorted-array-that-is-greater-than-the-target
        //https://www.geeksforgeeks.org/first-strictly-greater-element-in-a-sorted-array-in-java/
        int low = 0, high = intervals.getIntervalCount();
        int ans = -1;
        while (low <= high) {

            int mid = (low + high) / 2;
            if (intervals.getIntervalTime(mid) <= time) {
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
            ans = intervals.getIntervalCount() + 2;
        }
        return ans-1; // This finds the first interval that starts after the time so we want to shift by 1
    }

    @Override
    public int getEpochCount() {
        return epochCount;
    }

    @Override
    public double getEpochStartTime(int i) {
        return intervals.getIntervalTime(i);
    }

    @Override
    public double getEpochEndTime(int i) {
        return intervals.getIntervalTime(i) + intervals.getInterval(i);
    }

    @Override
    public double getEpochDuration(int i) {
        return intervals.getInterval(i);
    }

}
