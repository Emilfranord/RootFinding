
class HouseholderMod implements ItMe{ // Noor et al.: Modified Householder iterative method for nonlinear equations
	private NewtonRaphson nr;
	
	HouseholderMod(){
		nr = new NewtonRaphson();
	}
	public Double next(Double x, Func f){
		Double y = this.nr.next(x, f);
		Double fy = f.evaluate(y); 
		Func fp = f.differentiate();
		Func fpp = fp.differentiate();
		Double fpy = fp.evaluate(y);
		Double fppy = fpp.evaluate(y);
		
		return y - (fy / fpy) - ((fy*fy*fppy)/(2*fpy*fpy*fpy));
	}
	
	public String toString(){
		return "Modified Householder";
	}
}