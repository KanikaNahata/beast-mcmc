package dr.evomodel.coalescent.skybricks;

import dr.evolution.coalescent.*;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.evomodel.coalescent.CoalescentIntervalProvider;
import dr.evomodelxml.coalescent.skybrick.SkybrickLikelihoodParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.XMLParseException;
import org.apache.commons.math.util.FastMath;

import java.util.Arrays;
import java.util.List;

public class SkybrickLikelihood extends AbstractCoalescentLikelihood
        implements CoalescentIntervalProvider, Citable {


    private final Parameter popSizeParameter;
    private final EpochProvider epochProvider;
    private final IntervalList intervalsList;

    private final double[] sufficientStatistics;
    private final double[] numCoalEvents;

    protected int numEpochs;

    private final SkyBrickDemographicModel demographicModel;


    public SkybrickLikelihood(EpochProvider epochProvider, Parameter popSizes, IntervalList intervals, String demoType) throws XMLParseException {
        super(SkybrickLikelihoodParser.SKYBRICK,intervals);


        switch (demoType){
            case SkybrickLikelihoodParser.CONSTANT: demographicModel=new Constant();
                break;
            case SkybrickLikelihoodParser.EXPONENTIAL:demographicModel=new Exponential();
                break;
            default:
                throw new XMLParseException("Unrecognized demographic type");
        }

        this.intervalsList = intervals;
        this.popSizeParameter = popSizes;
        this.epochProvider = epochProvider;
        if(epochProvider instanceof Model){
            addModel((Model) epochProvider);
        }
        numEpochs = this.epochProvider.getEpochCount();
        this.numCoalEvents = new double[numEpochs];
        this.sufficientStatistics = new double[numEpochs];

        if(epochProvider instanceof FixedGridEpochProvider ){
            if(popSizeParameter.getDimension()!=numEpochs){
                throw new XMLParseException("Population size must be 1 greater than the number of grid points");
            }
        } else {
            if (demographicModel instanceof Exponential) {
                if (popSizeParameter.getDimension() != numEpochs + 1) {
                    throw new XMLParseException("Population size must be 1 greater than the number of epochs try: " + (numEpochs+1));
                }
            } else {
                if (popSizeParameter.getDimension() != numEpochs) {
                    throw new XMLParseException("Population size must be equal to the number of epochs, try: " + numEpochs);
                }
            }
        }

        addVariable(popSizes);
        calculateLogLikelihood();

    }



    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    @Override
    protected double calculateLogLikelihood() {
//        double coalescentLikelihood = piecewiseConstant ? calculateConstantCoalescentLogLikelihood() : calculateExpCoalescentLogLikelihood();
        // Matrix operations taken from block update sampler to calculate data likelihood and field prior
        setupSufficientStatistics();
        double currentLike = 0;
        double[] currentGamma = popSizeParameter.getParameterValues();

        for (int i = 0; i < numEpochs; i++) {
            currentLike += -numCoalEvents[i] * currentGamma[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
        }

        return currentLike;

    }

    protected void setupSufficientStatistics() {

        Arrays.fill(numCoalEvents, 0.0);
        Arrays.fill(sufficientStatistics, 0.0);
//        Arrays.fill(ploidySums, 0);
        //index of smallest grid point greater than at least one sampling/coalescent time in current tree
        int minEpoch;
        //index of greatest grid point less than at least one sampling/coalescent time in current tree
        int maxEpoch;
        int numLineages;
        int currentEpoch;
        int currentTimeIndex;

        double[] currentAndNextTime = new double[2];


        //time of last coalescent event in tree
        double lastCoalescentTime;

        currentTimeIndex = moveToNextTimeIndex(0, currentAndNextTime);

        //  numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
        numLineages = intervalsList.getLineageCount(currentTimeIndex);
        minEpoch = epochProvider.getEpoch(currentAndNextTime[0]);

        currentEpoch = minEpoch;

        lastCoalescentTime = currentAndNextTime[0] + intervalsList.getTotalDuration();

//            theLastTime = lastCoalescentTime;

        maxEpoch = epochProvider.getEpoch(lastCoalescentTime);

        SkyBrickDemographicModel demographic = getDemographicModel(currentEpoch);

        double currentEpochEndTime = epochProvider.getEpochEndTime(currentEpoch);
        double currentEpochStartTime = epochProvider.getEpochStartTime(currentEpoch);

        if (currentEpoch < maxEpoch) {

            while (currentAndNextTime[1] < currentEpochEndTime) {

                //check to see if interval ends with coalescent event
                //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                    numCoalEvents[currentEpoch] += demographic.getCoalescentEventWeight(currentEpochStartTime, currentAndNextTime[1]);
                }
                sufficientStatistics[currentEpoch] = sufficientStatistics[currentEpoch] + demographic.getScaledInterval(currentEpochStartTime, currentAndNextTime[0], currentAndNextTime[1]) * numLineages * (numLineages - 1) * 0.5;
                currentTimeIndex++;
                currentTimeIndex = moveToNextTimeIndex(currentTimeIndex, currentAndNextTime);
                numLineages = intervalsList.getLineageCount(currentTimeIndex);
            }
            sufficientStatistics[currentEpoch] = sufficientStatistics[currentEpoch] + demographic.getScaledInterval(currentEpochStartTime, currentAndNextTime[0], currentEpochEndTime) * numLineages * (numLineages - 1) * 0.5;
            currentEpoch++;


            while (currentEpoch < maxEpoch ) {

                //from likelihood of interval between first sampling time and minEpoch
                demographic = getDemographicModel(currentEpoch);
                currentEpochEndTime = epochProvider.getEpochEndTime(currentEpoch);
                currentEpochStartTime = epochProvider.getEpochStartTime(currentEpoch);

                if (currentAndNextTime[1] >= currentEpochEndTime) {
                    sufficientStatistics[currentEpoch] = sufficientStatistics[currentEpoch] + demographic.getScaledInterval(currentEpochStartTime, currentEpochStartTime, currentEpochEndTime) * numLineages * (numLineages - 1) * 0.5;
                    currentEpoch++;

                } else {
                    sufficientStatistics[currentEpoch] = sufficientStatistics[currentEpoch] + demographic.getScaledInterval(currentEpochStartTime, currentEpochStartTime, currentAndNextTime[1]) * numLineages * (numLineages - 1) * 0.5;
                    //check to see if interval ends with coalescent event
                    //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                    if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                        numCoalEvents[currentEpoch] += demographic.getCoalescentEventWeight(currentEpochStartTime, currentAndNextTime[1]);
                    }
                    currentTimeIndex++;
                    currentTimeIndex = moveToNextTimeIndex(currentTimeIndex, currentAndNextTime);
                    // numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
                    numLineages = intervalsList.getLineageCount(currentTimeIndex);

                    while (currentAndNextTime[1] < currentEpochEndTime && currentTimeIndex < intervalsList.getIntervalCount()) {
                        //check to see if interval ends with coalescent event
                        //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                        if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                            numCoalEvents[currentEpoch] += demographic.getCoalescentEventWeight(currentEpochStartTime, currentAndNextTime[1]);
                        }
                        sufficientStatistics[currentEpoch] = sufficientStatistics[currentEpoch] + demographic.getScaledInterval(currentEpochStartTime, currentAndNextTime[0], currentAndNextTime[1]) * numLineages * (numLineages - 1) * 0.5;
                        currentTimeIndex++;
                        currentTimeIndex = moveToNextTimeIndex(currentTimeIndex, currentAndNextTime);
                        numLineages = intervalsList.getLineageCount(currentTimeIndex);
                    }

                    sufficientStatistics[currentEpoch] = sufficientStatistics[currentEpoch] + demographic.getScaledInterval(currentEpochStartTime, currentAndNextTime[0], currentEpochEndTime) * numLineages * (numLineages - 1) * 0.5;
                    currentEpoch++;

                }

            }
            demographic = getDemographicModel(currentEpoch);
            currentEpochStartTime = epochProvider.getEpochStartTime(currentEpoch);
            sufficientStatistics[currentEpoch] = sufficientStatistics[currentEpoch] + demographic.getScaledInterval(currentEpochStartTime, currentEpochStartTime,currentAndNextTime[1]) * numLineages * (numLineages - 1) * 0.5;
            if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                numCoalEvents[currentEpoch] += demographic.getCoalescentEventWeight(currentEpochStartTime, currentAndNextTime[1]);
            }

            currentTimeIndex++;
            while ((currentTimeIndex) < intervalsList.getIntervalCount()) {

                currentTimeIndex = moveToNextTimeIndex( currentTimeIndex, currentAndNextTime);

                //numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
                numLineages = intervalsList.getLineageCount(currentTimeIndex);
                //check to see if interval is coalescent interval or sampling interval

                //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                    numCoalEvents[currentEpoch] += demographic.getCoalescentEventWeight(currentEpochStartTime, currentAndNextTime[1]);
                }
                sufficientStatistics[currentEpoch] = sufficientStatistics[currentEpoch] + demographic.getScaledInterval(currentEpochStartTime, currentAndNextTime[0], currentAndNextTime[1]) * numLineages * (numLineages - 1) * 0.5;
                currentTimeIndex++;

            }

        }else{
            // only one epoch
            while (currentTimeIndex < intervalsList.getIntervalCount()) {
                //check to see if interval ends with coalescent event
                //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                    numCoalEvents[currentEpoch] += demographic.getCoalescentEventWeight(currentEpochStartTime, currentAndNextTime[1]);
                }
                sufficientStatistics[currentEpoch] = sufficientStatistics[currentEpoch] + demographic.getScaledInterval(currentEpochStartTime, currentAndNextTime[0], currentAndNextTime[1]) * numLineages * (numLineages - 1) * 0.5;
                currentTimeIndex++;
                if ((currentTimeIndex) < intervalsList.getIntervalCount()) {
                    currentTimeIndex = moveToNextTimeIndex( currentTimeIndex, currentAndNextTime);

                    // numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
                    numLineages = intervalsList.getLineageCount(currentTimeIndex);
                }
            }

        }
    }



    private int moveToNextTimeIndex( int lastTimeIndex, double[] times) {
        int currentTimeIndex = lastTimeIndex;
        double currentTime = intervalsList.getIntervalTime(currentTimeIndex);
        double nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
        while (nextTime <= currentTime && currentTimeIndex + 2 < intervalsList.getIntervalCount()) {
            currentTimeIndex++;
            currentTime = intervalsList.getIntervalTime(currentTimeIndex);
            nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
        }
        times[0] = currentTime;
        times[1] = nextTime;
        return currentTimeIndex;
    }

    private SkyBrickDemographicModel getDemographicModel(int epoch){

        if(!(demographicModel instanceof Constant)) {

            SkyBrickDemographicModel demographicModel;
            if (epoch == epochProvider.getEpochCount() - 1 && epochProvider instanceof FixedGridEpochProvider) { // last epoch is always constant
                demographicModel = new Constant();
                return demographicModel;
            }
            // only need to scale for nonConstant population sizes
            this.demographicModel.setup(popSizeParameter.getParameterValue(epoch), popSizeParameter.getParameterValue(epoch + 1), epochProvider.getEpochDuration(epoch));
        }

        return this.demographicModel;

    }


    @Override
    public int getNumberOfCoalescentEvents() {
        return intervalsList.getIntervalCount() - intervalsList.getSampleCount() + 1;
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        throw new RuntimeException("What should this be?");
    }


    //Citable implementation
    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Skygrowth coalescent";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                new Author[]{
                        new Author("EM", "Volz"),
                        new Author("X", "Didelot"),
                },
                "Modeling the Growth and Decline of Pathogen Effective Population Size Provides Insight into Epidemic Dynamics and Drivers of Antimicrobial Resistance",
                2018,
                "Systematic Biology",
                67, 719, 728,
                "10.1093/sysbio/syy007"
        ));
    }

    /**
     * @return the units for this object.
     */
    @Override
    public Type getUnits() {
        return null;
    }

    /**
     * Sets the units for this object.
     *
     * @param units to use
     */
    @Override
    public void setUnits(Type units) {

    }


    /***
     * SkyBrick Demographic is Convert the timings and coalescent counts so that the sufficient statistics can be calculated as if
     * we were always dealing with a constant population size.
     */

    public interface SkyBrickDemographicModel {

        double getScaledInterval(double epochStateTime, double t0, double t1);

        double getCoalescentEventWeight(double epochStateTime, double t);


        void setup(double logN0, double logN1, double time);
    }

    class Constant implements SkyBrickDemographicModel{
            Constant(){}
            public void setup(double logN0, double logN1, double time){};
            public double  getScaledInterval(double epochStateTime, double t0, double t1){
                return t1-t0;
            };
            public double getCoalescentEventWeight(double epochStateTime,double t){
                return 1;
            };
    }

    class Exponential implements SkyBrickDemographicModel{
            Exponential(){}
            public void setup(double logN0, double logN1, double time){
                this.logN0 = logN0;

                this.r = (logN0 - logN1) / time;
            }
            private double getScaledTime(double t){
                return FastMath.exp(-r*t);
            }

            @Override
            public double getScaledInterval(double epochStateTime, double t0, double t1) {
                if(r==0){
                    return t1-t0;
                }
                // flipped to account for change of sign in exponential Log likelihood calculation compared to constant
                return (getScaledTime(t0 - epochStateTime) - getScaledTime(t1 - epochStateTime))/r ;
            }

            public double getCoalescentEventWeight(double epochStateTime,double t){
                if(logN0==0){
                    return 1; //doesn't matter since will be made 0 in calculate likelihoood
                }
                return (logN0 - (r * (t-epochStateTime)))/logN0 ;
            }
            private double logN0;
            private double r;
        }


}
