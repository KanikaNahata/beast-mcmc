package dr.evomodel.coalescent;

import dr.evolution.coalescent.*;
import dr.evomodelxml.coalescent.GMRFSkygrowthLikelihoodParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;
import dr.math.Binomial;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;
import org.apache.commons.math.util.FastMath;

import java.util.Arrays;
import java.util.List;

public class GMRFSkygrowthLikelihood extends AbstractCoalescentLikelihood
        implements CoalescentIntervalProvider, Citable {
    private final double[] gridPoints;
    private final int numGridPoints;


    private final Parameter popSizeParameter;
    private final Parameter precisionParameter;
    private final IntervalList intervalsList;

    private final ExponentialSkyGrowth demographicFunction = new ExponentialSkyGrowth();
    private final double[] growthRates;

    private final double[] sufficientStatistics;
    private final int[] numCoalEvents;
    protected SymmTridiagMatrix weightMatrix;


    private final GrowthRateStatistic grs;

    protected int fieldLength;


    public GMRFSkygrowthLikelihood(Parameter gridPoints, Parameter precision, Parameter popSizes, IntervalList intervals) {
        super(GMRFSkygrowthLikelihoodParser.SKYGROWTH,intervals);

        this.precisionParameter = precision;
        this.intervalsList = intervals;
        this.popSizeParameter = popSizes;
//        this.piecewiseConstant = piecewiseConstant;
//        if (popSizes.getDimension() != gridPoints.getDimension()) {

        if (popSizes.getDimension() - 1 != gridPoints.getDimension()) {
            throw new IllegalArgumentException("population sizes must be one greater than the number of grid points.");
        }

        this.gridPoints = gridPoints.getParameterValues();
        this.growthRates = new double[this.gridPoints.length];
        this.numGridPoints = this.gridPoints.length;

        fieldLength = this.gridPoints.length;

        this.numCoalEvents = new int[this.gridPoints.length];
        this.sufficientStatistics = new double[this.gridPoints.length];
        this.grs = new GrowthRateStatistic();

        setupGMRFWeights();
        addVariable(popSizes);
        addVariable(precision);

//        addModel((Model) intervalsList);
        calculateLogLikelihood();

    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************
    //TODO implement in abstract class?
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }


    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    @Override
    protected double calculateLogLikelihood() {
        double priorLikelihood = calculatePiorLogLikelihood();
//        double coalescentLikelihood = piecewiseConstant ? calculateConstantCoalescentLogLikelihood() : calculateExpCoalescentLogLikelihood();
        double coalescentLikelihood =  calculateExpCoalescentLogLikelihood();
        return coalescentLikelihood + priorLikelihood;
    }

    private DenseVector getCurrentRates() {
        double[] populationSizes = popSizeParameter.getParameterValues();

        //growth rates are defined forward in time but are indexed reverse in time
        for (int i = 0; i < growthRates.length; i++) {
            double h;
            if (i == 0) {
                h = gridPoints[i];
            } else {
                h = gridPoints[i] - gridPoints[i - 1];
            }
            growthRates[i] = ((populationSizes[i]) - (populationSizes[i + 1])) / h;
        }
        return new DenseVector(growthRates);
    }
    public SymmTridiagMatrix getScaledWeightMatrix(double precision) {
        SymmTridiagMatrix a = weightMatrix.copy();
        for (int i = 0; i < a.numRows() - 1; i++) {
            a.set(i, i, a.get(i, i) * precision);
            a.set(i + 1, i, a.get(i + 1, i) * precision);
        }
        a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
        return a;
    }
    protected void setupGMRFWeights() {

        //setupSufficientStatistics();

        //Set up the weight Matrix
        double[] offdiag = new double[fieldLength - 1];
        double[] diag = new double[fieldLength];

        //    private double theLastTime;
        double diagonalValue = 2;
        //First set up the offdiagonal entries;

        for (int i = 0; i < fieldLength - 1; i++) {
            offdiag[i] = -1;
        }

        //Then set up the diagonal entries;
        for (int i = 1; i < fieldLength - 1; i++) {
            //	diag[i] = -(offdiag[i] + offdiag[i - 1]);
            diag[i] = diagonalValue;
        }
        //Take care of the endpoints
        //diag[0] = -offdiag[0];
        //diag[fieldLength - 1] = -offdiag[fieldLength - 2];
        diag[0] = diagonalValue - 1.0;
        diag[fieldLength - 1] = diagonalValue - 1.0;


        weightMatrix = new SymmTridiagMatrix(diag, offdiag);

    }

    private double calculatePiorLogLikelihood() {

        DenseVector diagonal1 = new DenseVector(fieldLength);
        DenseVector currentRates = getCurrentRates();

        double currentLike = 0.0;

//        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0));
        currentQ.mult(currentRates, diagonal1);

        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentRates.dot(diagonal1);
//        if (lambdaParameter.getParameterValue(0) == 1) {
            currentLike -= (fieldLength - 1) / 2.0 * GMRFSkyrideLikelihood.LOG_TWO_TIMES_PI;
//        } else {
//            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
//        }
        return currentLike;



    }


    private double calculateExpCoalescentLogLikelihood() {
//index of smallest grid point greater than at least one sampling/coalescent time in current tree
        int minGridIndex;
        //index of greatest grid point less than at least one sampling/coalescent time in current tree
        int maxGridIndex;

        int numLineages;

        int currentGridIndex;
        int currentTimeIndex;

        double currentTime;
        double nextTime;
        double lastCoalescentTime;
//        boolean beyondTheWall = false;
        //numCoalEvents = new double[fieldLength];
        //sufficientStatistics = new double[fieldLength];

        double logL = 0;
        //from likelihood of interval between first sampling time and gridPoints[minGridIndex]
        currentTimeIndex = 0;
        currentTime = intervalsList.getIntervalTime(currentTimeIndex);
        nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
        while (nextTime <= currentTime) {
            currentTimeIndex++;
            currentTime = intervalsList.getIntervalTime(currentTimeIndex);
            nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
        }
        // need to reduce currentTimeIndex on getLineages

        //  numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
        numLineages = intervalsList.getLineageCount(currentTimeIndex);
        minGridIndex = 0;
        while (minGridIndex < numGridPoints - 1 && gridPoints[minGridIndex] <= currentTime) { // MAS: Unclear about need for -1
            minGridIndex++;
        }
        currentGridIndex = minGridIndex;

        lastCoalescentTime = currentTime + intervalsList.getTotalDuration();
//        maxGridIndex = numGridPoints - 1;
        maxGridIndex = numGridPoints - 2; //the last gridpoint is used to set the last rate but the interval goes back to the beginning of time
        while ((maxGridIndex >= 0) && (gridPoints[maxGridIndex] >= lastCoalescentTime)) {
            maxGridIndex = maxGridIndex - 1;
        }

        if (maxGridIndex >= 0 && minGridIndex < numGridPoints) {

            while (nextTime < gridPoints[currentGridIndex]) {

                this.demographicFunction.setup(popSizeParameter.getParameterValue(currentGridIndex),
                        popSizeParameter.getParameterValue(currentGridIndex + 1),
                        gridPoints[currentGridIndex]);
                //check to see if interval ends with coalescent event
                //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                    numCoalEvents[currentGridIndex]++;
                }
//                logL += calculateIntervalLikelihood(nextTime - currentTime, currentTime, numLineages, intervalsList.getIntervalType(currentTimeIndex));
                logL += calculateIntervalLikelihood(nextTime - currentTime, currentTime, numLineages, intervalsList.getIntervalType(currentTimeIndex));
//                sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - currentTime) * numLineages * (numLineages - 1) * 0.5;
                currentTime = nextTime;
                currentTimeIndex++;
                nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);

                while (nextTime <= currentTime) {
                    currentTimeIndex++;
                    currentTime = intervalsList.getIntervalTime(currentTimeIndex);
                    nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
                }

                //numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
                numLineages = intervalsList.getLineageCount(currentTimeIndex);

            }
            logL += calculateIntervalLikelihood(gridPoints[currentGridIndex] - currentTime, currentTime, numLineages, IntervalType.SAMPLE);

//            sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - currentTime) * numLineages * (numLineages - 1) * 0.5;

            currentGridIndex++;


            //from likelihood of intervals between gridPoints[minGridIndex] and gridPoints[maxGridIndex]

            this.demographicFunction.setup(popSizeParameter.getParameterValue(currentGridIndex),
                    popSizeParameter.getParameterValue(currentGridIndex + 1),
                    gridPoints[currentGridIndex] - gridPoints[currentGridIndex - 1]);

            while (currentGridIndex <= maxGridIndex) {
                if (nextTime >= gridPoints[currentGridIndex]) {

                    logL += calculateIntervalLikelihood((gridPoints[currentGridIndex] - gridPoints[currentGridIndex - 1]), 0, numLineages, IntervalType.SAMPLE);

                    currentGridIndex++;

                    this.demographicFunction.setup(popSizeParameter.getParameterValue(currentGridIndex),
                            popSizeParameter.getParameterValue(currentGridIndex + 1),
                            gridPoints[currentGridIndex] - gridPoints[currentGridIndex - 1]);
                } else {
                    logL += calculateIntervalLikelihood((nextTime - gridPoints[currentGridIndex - 1]), 0, numLineages, IntervalType.COALESCENT);

                    if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                        numCoalEvents[currentGridIndex]++;
                    }
                    currentTime = nextTime;
                    currentTimeIndex++;
                    nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
                    while (nextTime <= currentTime) {
                        currentTimeIndex++;
                        currentTime = intervalsList.getIntervalTime(currentTimeIndex);
                        nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
                    }

                    numLineages = intervalsList.getLineageCount(currentTimeIndex);

                    while (nextTime < gridPoints[currentGridIndex]) {
                        //check to see if interval is coalescent interval or sampling interval
                        if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {

                            numCoalEvents[currentGridIndex]++;
                        }

                        logL += calculateIntervalLikelihood(nextTime - currentTime, currentTime - gridPoints[currentGridIndex - 1], numLineages, intervalsList.getIntervalType(currentTimeIndex));

                        currentTime = nextTime;
                        currentTimeIndex++;
                        nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
                        while (nextTime <= currentTime) {
                            currentTimeIndex++;
                            currentTime = intervalsList.getIntervalTime(currentTimeIndex);
                            nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
                        }

                        numLineages = intervalsList.getLineageCount(currentTimeIndex);

                    }

                    logL += calculateIntervalLikelihood(gridPoints[currentGridIndex] - currentTime, currentTime - gridPoints[currentGridIndex - 1], numLineages, IntervalType.SAMPLE);
                    currentGridIndex++;

                    this.demographicFunction.setup(popSizeParameter.getParameterValue(currentGridIndex),
                            popSizeParameter.getParameterValue(currentGridIndex + 1),
                            gridPoints[currentGridIndex] - gridPoints[currentGridIndex - 1]);
                }
            }

            //from likelihood of interval between gridPoints[maxGridIndex] and lastCoalescentTime

            logL += calculateIntervalLikelihood((nextTime - gridPoints[currentGridIndex - 1]), 0, numLineages, IntervalType.COALESCENT);

            //check to see if interval ends with coalescent event
            // if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
            if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {

                numCoalEvents[currentGridIndex]++;
            }

            currentTime = nextTime;
            currentTimeIndex++;

            while ((currentTimeIndex) < intervalsList.getIntervalCount()) {

                nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
                while (nextTime <= currentTime) {
                    currentTimeIndex++;
                    currentTime = intervalsList.getIntervalTime(currentTimeIndex);
                    nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
                }

                numLineages = intervalsList.getLineageCount(currentTimeIndex);


                //check to see if interval is coalescent interval or sampling interval

                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
                    numCoalEvents[currentGridIndex]++;
                }

                logL += calculateIntervalLikelihood(nextTime - currentTime, currentTime - gridPoints[currentGridIndex - 1], numLineages, intervalsList.getIntervalType(currentTimeIndex));

                currentTime = nextTime;
                currentTimeIndex++;

            }
        }
        return logL;
    }


//    private double calculateConstantCoalescentLogLikelihood() {
//        //TODO only do this on tree move
//        setupSufficientStatistics();
//        // Matrix operations taken from block update sampler to calculate data likelihood and field prior
//        double currentLike = 0;
//        double[] currentGamma = popSizeParameter.getParameterValues();
//        for (int i = 0; i < numGridPoints; i++) {
//            currentLike += -numCoalEvents[i] * currentGamma[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
//        }
//        return currentLike;
//    }

//    protected void setupSufficientStatistics() {
//        //index of smallest grid point greater than at least one sampling/coalescent time in current tree
//        int minGridIndex;
//        //index of greatest grid point less than at least one sampling/coalescent time in current tree
//        int maxGridIndex;
//
//        int numLineages;
//
//        int currentGridIndex;
//        int currentTimeIndex;
//
//        double currentTime;
//        double nextTime;
//        double lastCoalescentTime;
//        //numCoalEvents = new double[fieldLength];
//        //sufficientStatistics = new double[fieldLength];
//
//        Arrays.fill(numCoalEvents, 0);
//        Arrays.fill(sufficientStatistics, 0);
//        //from likelihood of interval between first sampling time and gridPoints[minGridIndex]
//        currentTimeIndex = 0;
//        currentTime = intervalsList.getIntervalTime(currentTimeIndex);
//        nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
//        while (nextTime <= currentTime) {
//            currentTimeIndex++;
//            currentTime = intervalsList.getIntervalTime(currentTimeIndex);
//            nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
//        }
//        // need to reduce currentTimeIndex on getLineages
//
//        //  numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//        numLineages = intervalsList.getLineageCount(currentTimeIndex);
//        minGridIndex = 0;
//        while (minGridIndex < numGridPoints - 1 && gridPoints[minGridIndex] <= currentTime) { // MAS: Unclear about need for -1
//            minGridIndex++;
//        }
//        currentGridIndex = minGridIndex;
//
//        lastCoalescentTime = currentTime + intervalsList.getTotalDuration();
//        maxGridIndex = numGridPoints - 1;
//        while ((maxGridIndex >= 0) && (gridPoints[maxGridIndex] >= lastCoalescentTime)) {
//            maxGridIndex = maxGridIndex - 1;
//        }
//        //`if the last coalscent happens after the last grid point
//        // reject the state. This is a hard cut off
//        if (maxGridIndex == numGridPoints - 1) {
//            sufficientStatistics[0] = Double.POSITIVE_INFINITY;
//            return;
//        }
//
//        if (maxGridIndex >= 0 && minGridIndex < numGridPoints) {
//            while (nextTime < gridPoints[currentGridIndex]) {
//
//                //check to see if interval ends with coalescent event
//                //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
//                    numCoalEvents[currentGridIndex]++;
//                }
//                sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - currentTime) * numLineages * (numLineages - 1) * 0.5;
//                currentTime = nextTime;
//                currentTimeIndex++;
//                nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
//
//                while (nextTime <= currentTime) {
//                    currentTimeIndex++;
//                    currentTime = intervalsList.getIntervalTime(currentTimeIndex);
//                    nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
//                }
//
//                //numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//                numLineages = intervalsList.getLineageCount(currentTimeIndex);
//
//            }
//
//            sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - currentTime) * numLineages * (numLineages - 1) * 0.5;
//
//            currentGridIndex++;
//
//
//            //from likelihood of intervals between gridPoints[minGridIndex] and gridPoints[maxGridIndex]
//
//            while (currentGridIndex <= maxGridIndex) {
//                if (nextTime >= gridPoints[currentGridIndex]) {
//                    sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - gridPoints[currentGridIndex - 1]) * numLineages * (numLineages - 1) * 0.5;
//                    currentGridIndex++;
//                } else {
//
//                    sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - gridPoints[currentGridIndex - 1]) * numLineages * (numLineages - 1) * 0.5;
//
//                    //check to see if interval ends with coalescent event
//                    //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//                    if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
//                        numCoalEvents[currentGridIndex]++;
//                    }
//                    currentTime = nextTime;
//                    currentTimeIndex++;
//                    nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
//                    while (nextTime <= currentTime) {
//                        currentTimeIndex++;
//                        currentTime = intervalsList.getIntervalTime(currentTimeIndex);
//                        nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
//                    }
//
//                    // numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//                    numLineages = intervalsList.getLineageCount(currentTimeIndex);
//
//
//                    while (nextTime < gridPoints[currentGridIndex]) {
//                        //check to see if interval is coalescent interval or sampling interval
//                        //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//                        if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
//
//                            numCoalEvents[currentGridIndex]++;
//                        }
//                        sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - currentTime) * numLineages * (numLineages - 1) * 0.5;
//
//                        currentTime = nextTime;
//                        currentTimeIndex++;
//                        nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
//                        while (nextTime <= currentTime) {
//                            currentTimeIndex++;
//                            currentTime = intervalsList.getIntervalTime(currentTimeIndex);
//                            nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
//                        }
//
//                        //numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//                        numLineages = intervalsList.getLineageCount(currentTimeIndex);
//
//
//                    }
//                    sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (gridPoints[currentGridIndex] - currentTime) * numLineages * (numLineages - 1) * 0.5;
//
//                    currentGridIndex++;
//                }
//            }
//
//            //from likelihood of interval between gridPoints[maxGridIndex] and lastCoalescentTime
//
//            sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - gridPoints[currentGridIndex - 1]) * numLineages * (numLineages - 1) * 0.5;
//
//            //check to see if interval ends with coalescent event
//            // if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//            if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
//                numCoalEvents[currentGridIndex]++;
//            }
//
//            currentTime = nextTime;
//            currentTimeIndex++;
//
//            while ((currentTimeIndex) < intervalsList.getIntervalCount()) {
//                // currentTime = nextTime;
//                // currentTimeIndex++;
//
//                nextTime = intervalsList.getIntervalTime(currentTimeIndex+1);
//                while (nextTime <= currentTime) {
//                    currentTimeIndex++;
//                    currentTime = intervalsList.getIntervalTime(currentTimeIndex);
//                    nextTime = intervalsList.getIntervalTime(currentTimeIndex + 1);
//                }
//
//                //numLineages = intervalsList.get(i).getLineageCount(currentTimeIndex + 1);
//                numLineages = intervalsList.getLineageCount(currentTimeIndex);
//
//
//                //check to see if interval is coalescent interval or sampling interval
//
//
//                //if (intervalsList.get(i).getCoalescentEvents(currentTimeIndex + 1) > 0) {
//                if (intervalsList.getCoalescentEvents(currentTimeIndex) > 0) {
//                    numCoalEvents[currentGridIndex]++;
//                }
//                sufficientStatistics[currentGridIndex] = sufficientStatistics[currentGridIndex] + (nextTime - currentTime) * numLineages * (numLineages - 1) * 0.5;
//                currentTime = nextTime;
//                currentTimeIndex++;
//
//            }
//        }
//    }

    private double calculateIntervalLikelihood(
            double width, double timeOfPrevCoal, int lineageCount,
            IntervalType type) {
        final double timeOfThisCoal = width + timeOfPrevCoal;

        final double intervalArea = demographicFunction.getIntegral(timeOfPrevCoal, timeOfThisCoal);
        final double kchoose2 = Binomial.choose2(lineageCount);
        double like = -kchoose2 * intervalArea;

        switch (type) {
            case COALESCENT:
                final double demographic = demographicFunction.getLogDemographic(timeOfThisCoal);
                like += -demographic;

                break;
            case SAMPLE:
            case NOTHING:
                break;
        }

        return like;
    }

    @Override
    public int getNumberOfCoalescentEvents() {
        return intervalsList.getIntervalCount() - intervalsList.getSampleCount() + 1;
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        throw new RuntimeException("What should this be?");
    }

    public GrowthRateStatistic getGrowthRateStatistic() {
        return this.grs;
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


    public class GrowthRateStatistic extends Statistic.Abstract {

        /**
         * @return the number of dimensions that this statistic has.
         */
        @Override
        public int getDimension() {
            return growthRates.length;
        }

        /**
         * @param dim the dimension to return value of
         * @return the statistic's scalar value in the given dimension
         */
        @Override
        public double getStatisticValue(int dim) {
            return growthRates[dim];
        }
    }
    // A private version of ExponentialBSPGrowth that uses FastMath

    /**
     * This class models an exponentially growing (or shrinking) population
     * (Parameters: N0=present-day population size; r=growth rate).
     * This model is nested with the constant-population size model (r=0).
     *
     * @author Alexei Drummond
     * @author Andrew Rambaut
     * @version $Id: ExponentialGrowth.java,v 1.7 2005/05/24 20:25:56 rambaut Exp $
     */
    private class ExponentialSkyGrowth extends DemographicFunction.Abstract {

        /**
         * Construct demographic model with default settings
         */
        public ExponentialSkyGrowth() {
            super(Type.YEARS);
        }

        public void setup(double logN0, double logN1, double time) {
            this.N0 = FastMath.exp(logN0);
            this.r = (logN0 - logN1) / time;
        }

        // Implementation of abstract methods

        public double getDemographic(double t) {

            if (r == 0) {
                return N0;
            } else {
                return N0 * FastMath.exp(-t * r);
            }
        }

        public double getLogDemographic(double t) {
            if (r == 0) {
                return FastMath.log(N0);
            } else {
                return FastMath.log(N0) - (t * r);
            }
        }

        /**
         * Calculates the integral 1/N(x) dx between start and finish.
         */
        @Override
        public double getIntegral(double start, double finish) {
//        double integral1 = getNumericalIntegral(start, finish);

            double integral;
            if (r == 0.0) {
                integral = (finish - start) / N0;
            } else {
                integral = (FastMath.exp(finish * r) - FastMath.exp(start * r)) / N0 / r;
            }
            return integral;
        }

        /**
         * @return the number of arguments for this function.
         */
        public int getNumArguments() {
            return 0;
        }

        /**
         * @return the name of the n'th argument of this function.
         */
        public String getArgumentName(int n) {
            return null;
        }

        /**
         * @return the value of the n'th argument of this function.
         */
        public double getArgument(int n) {
            return 0;
        }

        /**
         * Sets the value of the nth argument of this function.
         */
        public void setArgument(int n, double value) {
        }

        /**
         * @return the lower bound of the nth argument of this function.
         */
        public double getLowerBound(int n) {
            return 0;
        }

        /**
         * Returns the upper bound of the nth argument of this function.
         */
        public double getUpperBound(int n) {
            return 0;
        }

        /**
         * Returns a copy of this function.
         */
        public DemographicFunction getCopy() {
            return null;
        }

        public double getIntensity(double t) {
            throw new RuntimeException("not implemented");
        }

        public double getInverseIntensity(double x) {
            throw new RuntimeException("not implemented");
        }


        private double r, N0;
    }

}
