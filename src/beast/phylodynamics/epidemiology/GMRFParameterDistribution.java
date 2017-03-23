package beast.phylodynamics.epidemiology;


import java.io.PrintStream;
import java.util.List;
import java.util.Random;


import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;



/**
 * Initial version Ported from Beast 1.7 ExponentialMarkovModel
 * Modified by David Rasmussen for a general prior on GMRF parameter
 */
@Description("A class that produces a distribution chaining values in a GMRF. ")
public class GMRFParameterDistribution extends Distribution {

    public Input<RealParameter> parameterInput = new Input<RealParameter>("parameter",
            "effective migration over time.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> groupSizesInput = new Input<RealParameter>("groupSizes",
            "group sizes (used to compute regime end point times).",Input.Validate.REQUIRED);
    
    public Input<RealParameter> gmrfPrecisionInput = new Input<RealParameter>("gmrfPrecision",
            "precision parameter for GMRF.",Input.Validate.REQUIRED);
	
    // **************************************************************
    // Private instance variables
    // **************************************************************
    //private RealParameter chainParameter = null;
    private double tau;

    @Override
    public void initAndValidate() throws Exception {
        
    	//chainParameter = parameterInput.get();;
    	tau = gmrfPrecisionInput.get().getValue();

    }


    /**
     * Get the log likelihood.
     *
     * @return the log likelihood.
     */
    @Override
    public double calculateLogP() throws Exception {
    	
    	logP = 0.0;
    	
		int regimes = groupSizesInput.get().getDimension();
		tau = gmrfPrecisionInput.get().getValue();
		double sumOfSqrDevs = 0.0;
		for (int r = 0; r < (regimes-1); r++) {
			final double dev = parameterInput.get().getArrayValue(r+1) - parameterInput.get().getArrayValue(r);
			sumOfSqrDevs += dev * dev;
		}
		logP -= 0.5 * tau * sumOfSqrDevs; // already in log units since we did not exponentiate
		logP += Math.log(Math.pow(tau, 0.5*(regimes-2)));
    	
        return logP;
        
    }
    
    /**
     * Loggable interface implementation follows *
     */
    @Override
    public void init(final PrintStream out) throws Exception {
    	
    	out.print(getID() + "\t");
    	//out.print("prior (" + getID() + ")\t");
    	
    }

    @Override
    public void log(final int nSample, final PrintStream out) {
    	
    	out.print(getCurrentLogP() + "\t");
        
    }

    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }

//    private double getChainValue(int i) {
//        if (uselog) {
//            return Math.log(chainParameter.getValue(index(i)));
//        } else {
//            return chainParameter.getValue(index(i));
//        }
//    }
//
//    private int index(int i) {
//        if (reverse)
//            return chainParameter.getDimension() - i - 1;
//        else
//            return i;
//    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }
}

