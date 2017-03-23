package beast.phylodynamics.model;

import java.util.ArrayList;
import java.util.Collections;

import org.jblas.DoubleMatrix;

import beast.core.CalculationNode;
import beast.core.Description;

/**
 * @author Nicola Felix Mueller
 */
@Description("Saves the Timeseries outside of the Epidemiological class" +
			 "Only used to increase overview in the Epi class")
public class TimeSeries extends CalculationNode {

    public ArrayList<Double> integrationTimes;						// times at which euler integration is performed
    public ArrayList<DoubleMatrix> fSeries; 						// Transmission rates
    public ArrayList<DoubleMatrix> gSeries; 						// Any transition that does not need transmission
    public ArrayList<DoubleMatrix> ySeries;							// population sizes over time 

    public ArrayList<DoubleMatrix> mSeries;							// moving Fractions over time 

    public ArrayList<DoubleMatrix> cSeries;							// cumulative incidence over time
    
    public ArrayList<DoubleMatrix> effRSeries;						// effective reproductive number over time
    
    public ArrayList<DoubleMatrix> nSTraj;							// Trajectory of susceptibles
    public ArrayList<DoubleMatrix> nITraj;							// Trajectory of infected
    public ArrayList<DoubleMatrix> nRTraj;							// Trajectory of infected
	
	@Override
	public void initAndValidate() throws Exception {

	}
	
	/**
	 * set Back all the ArrayLists and DoubleMatrices
	 * @param states
	 */
	protected void setBack(int states){
		integrationTimes = new ArrayList<Double>();
		
		fSeries = new ArrayList<DoubleMatrix>();
		gSeries = new ArrayList<DoubleMatrix>();
		ySeries = new ArrayList<DoubleMatrix>();
		mSeries = new ArrayList<DoubleMatrix>();
		
		cSeries = new ArrayList<DoubleMatrix>();
		effRSeries = new ArrayList<DoubleMatrix>();
		
		nSTraj = new ArrayList<DoubleMatrix>();
		nITraj = new ArrayList<DoubleMatrix>();
		nRTraj = new ArrayList<DoubleMatrix>();
		
		currS = new DoubleMatrix(states);
		currI = new DoubleMatrix(states);
		currR = new DoubleMatrix(states);
		currN = new DoubleMatrix(states);
		printS = new DoubleMatrix(states);
		printI = new DoubleMatrix(states);
		
	}

	
	/*
	 * Add to the ArrayLists
	 */
	
	protected void addS(DoubleMatrix S){
		DoubleMatrix s = new DoubleMatrix();
		s.copy(S);
		nSTraj.add(s);
		s = null;
	}
	
	protected void addI(DoubleMatrix I){
		DoubleMatrix i = new DoubleMatrix();
		i.copy(I);
		nITraj.add(i);
		i = null;
	}
	
	protected void addR(DoubleMatrix R){
		DoubleMatrix r = new DoubleMatrix();
		r.copy(R);
		nRTraj.add(r);
		r = null;
	}
	
	protected void addF(DoubleMatrix F){
		DoubleMatrix f = new DoubleMatrix();
		f.copy(F);
		fSeries.add(f);
		f = null;
	}
	
	protected void addG(DoubleMatrix G){
		DoubleMatrix g = new DoubleMatrix();
		g.copy(G);
		gSeries.add(g);
		g = null;
	}
	
	protected void addY(DoubleMatrix Y){
		DoubleMatrix y = new DoubleMatrix();
		y.copy(Y);
		ySeries.add(y);
		y = null;
	}
	
	protected void addM(DoubleMatrix M){
		DoubleMatrix m = new DoubleMatrix();
		m.copy(M);
		mSeries.add(m);
		m = null;
	}
	
	protected void addIncidence(DoubleMatrix C){
		DoubleMatrix c = new DoubleMatrix();
		c.copy(C);
		cSeries.add(c);
		c = null;
	}
	
	protected void addEffR(DoubleMatrix EffR){
		DoubleMatrix effR = new DoubleMatrix();
		effR.copy(EffR);
		effRSeries.add(effR);
		effR = null;
	}
	

	/*
	 * Get Values From the ArrayLists
	 */
	protected double getT(int i){
		return integrationTimes.get(i);
	}
	
	protected DoubleMatrix getS(int i){
		return nSTraj.get(i);
	}
	
	protected DoubleMatrix getI(int i){
		return nITraj.get(i);
	}
	
	protected DoubleMatrix getR(int i){
		return nRTraj.get(i);
	}
	
	protected DoubleMatrix getF(int i){
		return fSeries.get(i);
	}
	
	protected DoubleMatrix getG(int i){
		return gSeries.get(i);
	}
	
	protected DoubleMatrix getY(int i){
		return ySeries.get(i);
	}
	
	protected DoubleMatrix getM(int i){
		return mSeries.get(i);
	}
	
	protected DoubleMatrix getC(int i){
		return cSeries.get(i);
	}
	
	protected DoubleMatrix getEffR(int i){
		return effRSeries.get(i);
	}
	
	protected int size(){
		return nSTraj.size();
	}
	
	
	
	/*
	 * Reverse the ArrayLists
	 */
	protected void reverse(){
		
		Collections.reverse(fSeries);
		Collections.reverse(gSeries);
		Collections.reverse(ySeries);
		Collections.reverse(mSeries);
		Collections.reverse(cSeries);
		Collections.reverse(effRSeries);
		Collections.reverse(nSTraj);
		Collections.reverse(nITraj);
		Collections.reverse(nRTraj);
		Collections.reverse(integrationTimes);
		
	}
	
    public DoubleMatrix currS;
    public DoubleMatrix currI;
    public DoubleMatrix currR;
    public DoubleMatrix currN;
    public DoubleMatrix printI;
    public DoubleMatrix printS;
	
	/*
	 * Set the Double Matrices
	 */
	protected void setS(DoubleMatrix S){
		currS.copy(S);
	}
	
	protected void setI(DoubleMatrix I){
		currI.copy(I);
	}
	
	protected void setI(int i, double init){
		currI.put(i, init);
	}
	
	protected void setS(int i, double init){
		currS.put(i, init);
	}
	
	protected void setR(DoubleMatrix R){
		currR.copy(R);
	}
	
	protected void setN(DoubleMatrix N){
		currN.copy(N);
	}
	
	protected void setPS(DoubleMatrix S){
		printI.copy(S);
	}
	
	/*
	 * get the Double Matrices
	 */
	protected DoubleMatrix getS(){
		return currS;
	}
	
	protected DoubleMatrix getI(){
		return currI;
	}
	
	protected DoubleMatrix getR(){
		return currR;
	}
	
	protected DoubleMatrix getN(){
		return currN;
	}
	
	protected DoubleMatrix getPS(){
		return printI;
	}
	
	protected DoubleMatrix getPI(){
		return printS;
	}
	


}
