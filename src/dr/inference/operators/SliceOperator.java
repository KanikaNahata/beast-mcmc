package dr.inference.operators;

import dr.inference.prior.Prior;
import dr.inference.model.*;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.util.Attribute;
import dr.math.distributions.NormalDistribution;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements a generic univariate slice sampler.
 *
 * See: RM Neal (2003) Slice Sampling, Annals of Statistics, 31, 705-767 (with discussion)
 *
 * @author Marc A. Suchard
 */
public class SliceOperator extends SimpleMetropolizedGibbsOperator {

    public SliceOperator(Variable<Double> variable) {
        this(new SliceInterval.SteppingOut(), variable);
    }

    public SliceOperator(SliceInterval sliceInterval, Variable<Double> variable) {
        this.sliceInterval = sliceInterval;

        if (variable.getSize() != 1) {
            throw new RuntimeException("Generic slice sampler is currently for univariate parameters only");
        }        
        this.variable = variable;
    }

    public double doOperation(Prior prior, Likelihood likelihood) throws OperatorFailedException {
        double newValue = sliceInterval.drawFromInterval(this, 1.0, width);
        variable.setValue(0, newValue);
        return 0;
    }

    public int getStepCount() {
        return 1;
    }

    public String getOperatorName() {
        return "genericSliceSampler";
    }

    public static void main(String[] arg) {

        // Define normal model
        Parameter mean = new Parameter.Default(1.0); // Starting value
        Variable<Double> stdev = new Variable.D(1.0, 1); // Fixed value
        ParametricDistributionModel densityModel = new NormalDistributionModel(mean, stdev);
        DistributionLikelihood likelihood = new DistributionLikelihood(densityModel);

        // Define prior
        DistributionLikelihood prior = new DistributionLikelihood(new NormalDistribution(0.0, 1.0)); // Hyper-priors
        prior.addData(mean);

        // Define data
        likelihood.addData(new Attribute.Default<double[]>("Data", new double[] {0.0, 1.0, 2.0}));

        List<Likelihood> list = new ArrayList<Likelihood>();
        list.add(likelihood);
        list.add(prior);
        CompoundLikelihood posterior = new CompoundLikelihood(0, list);
        SliceOperator sliceSampler = new SliceOperator(mean);

        final int length = 1000;
        double total = 0;

        for(int i = 0; i < length; i++) {
            try {
                sliceSampler.doOperation(null, posterior);
            } catch (OperatorFailedException e) {
                System.err.println(e.getMessage());
            }
            total += mean.getValue(0);
        }
        System.out.println("E(x) = "+(total/length));
    }

    private final SliceInterval sliceInterval;
    private final double width = 1.0;
    private final Variable<Double> variable;
}
