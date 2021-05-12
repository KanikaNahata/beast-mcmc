package dr.evomodel.coalescent.skybricks;

public interface EpochProvider {
    int getEpoch(double time);

    int getEpochCount();

    double getEpochStartTime(int i);

    double getEpochEndTime(int i);

    double getEpochDuration(int i);
}
