
class HalleyMod implements ItMe{ // Noor et al.: A new modified Halley method without second derivatives for nonlinear equation
	private NewtonRaphson nr;
	
	HalleyMod(){
		nr = new NewtonRaphson();
	}

	public Double next(Double x, Func f){
		Double fx  =  f.evaluate(x);
		Func fp = f.differentiate();
		Double fpx =  fp.evaluate(x);
		Double y   =  this.nr.next(x, f);
		Double fy  =  f.evaluate(y);
		Double fpy =  fp.evaluate(y);
		return y - ((2*fx*fy*fpy)/(2*fx*fpy*fpy - fpx*fpx*fy + fpx*fpy*fy));
	}
	public String toString(){
		return "Modified Halley";
	}
}