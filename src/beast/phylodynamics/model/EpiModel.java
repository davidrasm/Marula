package beast.phylodynamics.model;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math.MathException;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.ParametricDistribution;


/**
 * @author David Rasmussen
 */
@Description("Epidemic model" +
             "Tracks deterministic pop trajectory and F, G and Y matrices")
public class EpiModel extends CalculationNode implements Loggable {
	
    public Input<RealParameter> timeStepInput = new Input<RealParameter>("timeStep",
            "timeStep (defaults to 2).",Input.Validate.REQUIRED);    

    /*
     * Input specifying the origin, meaning the point at which to start with
     * simulating the S(E)IR trajectories as well as the f- g- and y-series 
     */
    public Input<RealParameter> originInput = new Input<RealParameter>("origin",
            "timepoint of origin.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> popNInput = new Input<RealParameter>("popN",
            "nodes in network.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> tauInput = new Input<RealParameter>("tau",
            "transmission rate.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> nuInput = new Input<RealParameter>("nu",
            "removal rate.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> noImmunityInput = new Input<RealParameter>("noImmunity",
            "immunity after infection.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> birthRateInput = new Input<RealParameter>("birthRate",
            "birth rate.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> deathRateInput = new Input<RealParameter>("deathRate",
            "death rate.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> gammaInput = new Input<RealParameter>("gamma",
            "age transition rate.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> betaInput = new Input<RealParameter>("beta",
            "beta.",Input.Validate.REQUIRED);
    
    public Input<BooleanParameter> betaIndicatorInput = new Input<BooleanParameter>("betaIndicator",
            "beta.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> endPointTimesInput = new Input<RealParameter>("endPointTimes",
            "end point times.",Input.Validate.REQUIRED);
    
    public Input<RealParameter> groupSizesInput = new Input<RealParameter>("groupSizes",
            "group sizes (used to compute regime end point times).",Input.Validate.REQUIRED);
    
    public Input<RealParameter> R0Input = new Input<RealParameter>("R0",
            "net reproductive rate R0.");
    
    public Input<IntegerParameter> statesInput = new Input<IntegerParameter>(
    		"infectedStates",
    		"total number of infected classes");
    
    /*
     * Boolean Inputs used to define which model to use. Would technically not 
     * be necessary but is nice to have
     */    
    public Input<BooleanParameter> isSIRInput = new Input<BooleanParameter>(
    		"isSIR",
    		"define whether this is an SIS or SIR model");
    
    public Input<BooleanParameter> useEndPointTimesInput = new Input<BooleanParameter>(
    		"useEndPointTimes",
    		"define whether to use input EndPointTimes or set by groupSizes input");
    
    public Input<Boolean> estimateR0Input = new Input<Boolean>(
    		"estimateR0",
    		"if estimating R0 and computing tau based on R0. Default is false", false);

    /*
     * Parameters that are needed independent of the population model used
     */
    public boolean dirty;
    public boolean reject = false;
    public boolean diagTrans = false;    
    
    /*
     * TimeSeries Object that stores all the F, G, Y and S I R ArraLists
     */
    public TimeSeries timeSeries;
    
    /*
     * ArrayList storing all the times between which euler integration is
     * performed
     */
    public ArrayList<Double> integrationTimes;
   
    public int states;
    
    // Parameter vectors
    protected DoubleMatrix nuVec;
    protected DoubleMatrix noImmunityVec;
    protected DoubleMatrix popNVec;
    protected DoubleMatrix birthRateVec;
    protected DoubleMatrix deathRateVec;
    protected DoubleMatrix gammaVec;
    protected DoubleMatrix betaMatrix;
    protected DoubleMatrix cumlIncidence;
    
    protected DoubleMatrix endPointTimes;
    protected int regimes;
   
    /*
     * Boolean that store the population model used for 
     * calculating the F G & Y matrices over time. At least
     * one of those has to be true 
     */
    
    public boolean SIR = false;
    public boolean useEndPointTimes = false;
    public boolean estimateR0 = false;
    
    private double origin;
    private double timeStep;
    	
	@Override
	public void initAndValidate() throws Exception {
		
		states = statesInput.get().getValue();
		timeStep = timeStepInput.get().getValue();
		
		if(originInput.get() != null) origin  = originInput.get().getValue();
		
		/*
		 * Get the population model used as Input
		 */
		if (isSIRInput.get() != null) SIR = isSIRInput.get().getValue();
		
		/*
		 * Set usage of GMRF model
		 */
		if (useEndPointTimesInput.get() != null) useEndPointTimes = useEndPointTimesInput.get().getValue();
		
		/*
		 * If R0 is estimated than tau will be computed from R0
		 */
		if (estimateR0Input.get() != null) estimateR0 = estimateR0Input.get();

		timeSeries = new TimeSeries();
	}
	
	public boolean update() throws MathException{
		
		timeSeries.setBack(states);	// empties all arrays in timeS
		
    	boolean reject = false;
    	
    	if (!useEndPointTimes) {
    		
    		setEndPointTimes();
    		
    	} else {
    	
    		// Check ordering of endPointTimes
    		regimes = endPointTimesInput.get().getDimension();
    		endPointTimes = new DoubleMatrix(regimes);
    		double previousTime = origin;
    		for (int r = 0; r < regimes; r++) {
    			double nextRegimeTime;
    			if (r == regimes-1) {
    				nextRegimeTime = 0.0;
    			} else {
    				nextRegimeTime = endPointTimesInput.get().getValue(r);
    			}
    			if (nextRegimeTime >= previousTime) {
    				System.out.println("Bad ordering");
    				reject = true;
    				return reject;
    			}
    			endPointTimes.put(r, nextRegimeTime);
    			previousTime = nextRegimeTime;	
    		}
    		
    	}
    	
//      if (!dirty) return reject;        
        
        final double endTime = 0.0;
        
        // Populate parameter vectors without secular trends
        nuVec = new DoubleMatrix(states);
        noImmunityVec = new DoubleMatrix(states);
        popNVec = new DoubleMatrix(states);
        birthRateVec = new DoubleMatrix(states);
        deathRateVec = new DoubleMatrix(states);
        gammaVec = new DoubleMatrix(states); // if we have migration other than through transmission
        for (int k = 0; k < states; k++) {
        	nuVec.put(k, nuInput.get().getArrayValue(k));
        	noImmunityVec.put(k, noImmunityInput.get().getArrayValue(k));
        	popNVec.put(k, popNInput.get().getArrayValue(k));
        	birthRateVec.put(k, birthRateInput.get().getArrayValue(k));
        	deathRateVec.put(k, deathRateInput.get().getArrayValue(k));
        	gammaVec.put(k, gammaInput.get().getArrayValue(k));
        }
        
        // Construct beta matrix
        DoubleMatrix beta = DoubleMatrix.zeros(states,states);
        int linearIndex = 0;
        int regime = 0;
        for (int k = 0; k < states; k++) {
        	for (int l = 0; l < states; l++) {
        		linearIndex = (k * states * regimes) + (l * regimes) + regime; // linear index in paramInput array
        		beta.put(k, l, betaInput.get().getArrayValue(linearIndex));
        	}
        }
        
        // Get initial state variables
        DoubleMatrix initI = DoubleMatrix.zeros(states,1);
        initI.put(0, 1.0); // set initial pop size in global pop to 1.0
        DoubleMatrix initS = popNVec.sub(initI);
        cumlIncidence = DoubleMatrix.zeros(states,states);
        
        // Set time series at t = 0
        timeSeries.addIncidence(cumlIncidence);
        timeSeries.setS(initS);
        timeSeries.setI(initI);
        timeSeries.setN(popNVec);
        
        // Push init values to lists
		timeSeries.addS(timeSeries.getS());
		timeSeries.addI(timeSeries.getI());
		timeSeries.addEffR(beta.diag().mul(timeSeries.getS().div(timeSeries.getN())).div(nuVec));
        
        // Set integrationTimes for Euler integration
        integrationTimes = new ArrayList<Double>();    	
		double time = origin;
		while (time > endTime){
			integrationTimes.add(time);
			time -= timeStep;
		}		
		integrationTimes.add(endTime);

		
		int t = 0;
		int dur = integrationTimes.size();
		
		// If one regime
		//do{
		//	reject = deltaStep(t, dur, timeStep);
		//	t++;
		//} while(integrationTimes.get(t-1)>0 && !reject);
		
		// If more than one regime
		double nextRegimeTime = endPointTimes.get(0);
		for (int r = 0; r < regimes; r++) {
			
			DoubleMatrix betaRegime = new DoubleMatrix();
			betaRegime.copy(beta);
	        for (int k = 0; k < states; k++) {
	        	for (int l = 0; l < states; l++) {
	        		linearIndex = (k * states * regimes) + (l * regimes) + r; // linear index in paramInput array
	        		if (betaIndicatorInput.get().getArrayValue(linearIndex) > 0) { 
	        			// If we are estimating this parameter, otherwise leave at base value
	        			betaRegime.put(k, l, betaInput.get().getArrayValue(linearIndex));
	        		}
	        	}
	        }
			betaMatrix = betaRegime;
			
			if (r == regimes-1) {
				nextRegimeTime = 0.0;
			} else {
				nextRegimeTime = endPointTimes.get(r);
			}

			do{
				reject = deltaStep(t, dur, timeStep, r);
				t++;
			} while(integrationTimes.get(t-1)>nextRegimeTime && !reject);

		}

		/*
		 * Reverse the array Lists in order to go from
		 * forward in time to backwards in time
		 */
		timeSeries.reverse();
		Collections.reverse(integrationTimes);

		return reject;		
	}
	
	/*
	 * Calculation of the exponential trajectories as well as updating of
	 * the initial time point of the first individual
	 */
	protected boolean deltaStep(int t, int dtTimeCount, double TimeStep, int regime){
		
		// Transmission
		DoubleMatrix sRowVec = timeSeries.getS().div(timeSeries.getN()).transpose();
	    DoubleMatrix transMatrix = betaMatrix.mul(timeSeries.getI().mmul(sRowVec)); // beta * I * S
	    
	    // Create copy to store transMatrix for coalescent rates
	    DoubleMatrix Transmission = new DoubleMatrix(); // transmission is from k --> l
	    Transmission.copy(transMatrix);
	    
	    // For state variables, set transmission from influxLocs -> external pop to zero
	    double transExternal = transMatrix.get(0,0);
	    DoubleMatrix newC = DoubleMatrix.zeros(states,1);
	    transMatrix.putColumn(0, newC); // to keep external pop in equilibrium
	    transMatrix.put(0,0, transExternal);
	    DoubleMatrix transVec = transMatrix.columnSums().transpose();  // sum is over rows
		
		// Removal
	    DoubleMatrix recVec = timeSeries.getI().mul(nuVec);
	    
	    // Births
	    //DoubleMatrix birthVec = birthRateVec; // as a number
	    DoubleMatrix birthVec = timeSeries.getN().mul(birthRateVec); // as a per capita rate
	    
	    // Deaths
	    DoubleMatrix deathSVec = timeSeries.getS().mul(deathRateVec);
	    DoubleMatrix deathIVec = timeSeries.getI().mul(deathRateVec);
	    
	    // Age-transitions (new way)
	    DoubleMatrix currS = timeSeries.getS();
	    DoubleMatrix migMatrixS = DoubleMatrix.zeros(states, states);;
	    DoubleMatrix migOutSVec = migMatrixS.rowSums();
	    DoubleMatrix migInSVec = migMatrixS.columnSums().transpose();
	    
	    DoubleMatrix currI = timeSeries.getI();
	    DoubleMatrix migMatrixI = DoubleMatrix.zeros(states, states);
	    DoubleMatrix migOutIVec = migMatrixI.rowSums();
	    DoubleMatrix migInIVec = migMatrixI.columnSums().transpose();
	    
	    // Change in state variables
	    DoubleMatrix dS = birthVec.add(noImmunityVec.mul(recVec)).sub(transVec).sub(deathSVec).sub(migOutSVec).add(migInSVec);
	    DoubleMatrix dI = transVec.sub(recVec).sub(deathIVec).sub(migOutIVec).add(migInIVec);
		
		newC = Transmission.getColumn(0).div(sRowVec); // to get beta * I
		Transmission.putColumn(0, newC);
		Transmission.put(0,0,transExternal);
		
		//DoubleMatrix Transmission = transMatrix; // transmission is from k --> l
		DoubleMatrix Migration = migMatrixI; // migration from k --> l
		
		// Have to add F, G and Y
		timeSeries.addF(Transmission);
		timeSeries.addG(Migration);
		timeSeries.addY(timeSeries.getI());
		
		if (t < (dtTimeCount - 1)) {
		    
	    	// Update state variables
	    	final double dt = integrationTimes.get(t) - integrationTimes.get(t+1);
	    		
	    	//DoubleMatrix newS = timeS.getS().add(dS.mul(dt));
	    	timeSeries.setS(timeSeries.getS().add(dS.mul(dt)));
	    	//DoubleMatrix newI = timeS.getI().add(dI.mul(dt));
	    	timeSeries.setI(timeSeries.getI().add(dI.mul(dt)));
	    	
	    	// Update pop densities ???
	    	//timeS.setN(timeS.getS().add(timeS.getI()));
	    	
			// Add cumulative incidence
			cumlIncidence = cumlIncidence.add(Transmission.mul(dt));
			timeSeries.addIncidence(cumlIncidence);
	    	
	    	if (timeSeries.getPS().min() < 0|| timeSeries.getPI().min() < 0){
	    		return true;
	    	}
           	timeSeries.addS(timeSeries.getS());
           	timeSeries.addI(timeSeries.getI());
           	timeSeries.addEffR(betaMatrix.diag().mul(timeSeries.getS().div(timeSeries.getN())).div(nuVec));
           	
           	// Consistancy checks
           	//double indvCount = timeS.getS().sum() + timeS.getI().sum();
           	//currI = timeS.getI();
           	//currS = timeS.getS();
           	//System.out.println(integrationTimes.get(t+1));
           	//System.out.println(currI);
           	//System.out.println(currS);
           	//System.out.println(indvCount);
           	
	    }
		
		return false;
	}
	
	public void setEndPointTimes() {
    	
		// Compute end pointTimes
    	regimes = groupSizesInput.get().getDimension();
    	endPointTimes = new DoubleMatrix(regimes);
		double previousTime = origin;
		for (int r = 0; r < regimes; r++) {
			double nextRegimeTime;
			if (r == regimes-1) {
				nextRegimeTime = 0.0;
			} else {
				nextRegimeTime = previousTime - groupSizesInput.get().getValue(r);
				if (nextRegimeTime < 0) {
					nextRegimeTime = 0.0;
					System.out.println("Sum of group size times exceeds maximum possible time");	
				}
			}
			endPointTimes.put(r, nextRegimeTime);
			previousTime = nextRegimeTime;
		}
	}
	
	protected double computeTauFromR(DoubleMatrix Nkl, DoubleMatrix Nk) {
		
		double tRate = 0;
		return tRate;
		
	}
	
	/*
	private void updateInitialTimePoint(){
		intitalInfected = new DoubleMatrix(states);
		for (int s = 0; s < states; s++)
			intitalInfected.put(s, intitalIntroductionInput.get().getValue(s));
	} */
	
	/*
	 * something else
	 */
	protected boolean delta(int i, int t, double s){
		return false;
	}
	
	public int getTrajLength(){
    	return timeSeries.size();
    }

    public DoubleMatrix getNumS(int t) {
    	return timeSeries.getS(t);
    }
    
	public DoubleMatrix getNumI(int t) {
		return timeSeries.getI(t);
    }
    
	public DoubleMatrix getNumR(int t) {
		return timeSeries.getR(t);
    }
	
	public double getNumSDeme(int t, int d) {
		return timeSeries.getS(t).get(d);
    }
    
	public double getNumIDeme(int t, int d) {
		return timeSeries.getI(t).get(d);
    }
    
	public double getNumRDeme(int t, int d) {
		return timeSeries.getR(t).get(d);
    }
	
	public double getTime(int t) {
		return integrationTimes.get(t);
    }
    
    /**
     * @param t time interval to return
     * @param k state to return Y for
    */
    public double getY(int t, int k) {
    	return timeSeries.getY(t).get(k);
    }
    
    public DoubleMatrix getY(int t) {
    	return timeSeries.getY(t);
    }
    
    /**
     * @param t time interval to return
     * @param k row in F
     * @param l column in F
    */
    public double getF(int t, int k, int l) {
    	return timeSeries.getF(t).get(k,l);
    }
    
    public DoubleMatrix getF(int t) {
    	return timeSeries.getF(t);
    }
    
    /**
     * @param t time interval to return
     * @param k row in G
     * @param l column in G
    */
    public double getG(int t, int k, int l) {
    	return timeSeries.getG(t).get(k,l);
    }
    
    public DoubleMatrix getG(int t) {
    	return timeSeries.getG(t);
    }
    
    /**
     * 
     * @param t
     * @param k
     * @param l
     * @return
     */
    public double getM(int t, int k, int l) {
    	return timeSeries.getM(t).get(k,l);
    }
    
    public DoubleMatrix getM(int t) {
    	return timeSeries.getM(t);
    }
    
    /**
     * @param t time interval to return
     * @param k row in F
     * @param l column in F
    */
    public double getC(int t, int k, int l) {
    	return timeSeries.getC(t).get(k,l);
    }
    
    public double getEffR(int t, int k) {
    	return timeSeries.getEffR(t).get(k);
    }
    
    /**
     * @return integrationTimes in reverse order (backwards in time)
     */
    public ArrayList<Double> getIntegrationTimes() {
    	return integrationTimes;
    }
    
    public double getOrigin(){
    	return origin;
    }
    
    public void setOrigin(double origin){
    	this.origin = origin;
    }
    
    /*
     * CalculationNode interface
     */

    @Override
    public boolean requiresRecalculation() {
        dirty = true;
        return true;
    }

    @Override
    public void restore() {
        dirty = true;
        super.restore();
    }
    
    public void init(PrintStream out) throws Exception {

        //out.print("meanK\t");
        //out.print("varK\t");
        //out.print("tau\t");

    }

    public void log(int nSample, PrintStream out) {
    	
        //out.format("%g\t", meanK);
    	//out.format("%g\t", varK);
    	//out.format("%g\t", tau);
    	//System.out.println(pk);
        
    	// was this:
        //out.format("%g\t", betaParameter.get().getValue()*n_S_Parameter.get().getValue()
                ///gammaParameter.get().getValue());

        //double tend = NStraj.size() * dt;
        //double delta = tend / (statesToLogInput.get() - 1);

    }

    public void close(PrintStream out) {
    }
    
	
	
}

