
class DecompositionIII implements ItMe{ //Chun, C.:Iterative methods improving Newton's method by the decomposition method 
	private NewtonRaphson nr;
	DecompositionIII(){
		nr = new NewtonRaphson();
	}
	public Double next(Double x, Func f){
		Double fx  = f.evaluate(x);
		Func fp = f.differentiate();
		Func fpp = fp.differentiate();
		Double fpx = fp.evaluate(x);
		Double y   = this.nr.next(x, f);
		Double fy  = f.evaluate(y);
		Double fpy = fp.evaluate(y);		
		Double fppy = fpp.evaluate(y);
		
		return x - (fx)/(fpx) - (3.0*fy)/(fpx) + (3.0*fy*fpy)/(fpx*fpx) - (0.5 * (fy * (fy*fppy+2.0*fpy*fpy)))/(fpy*fpy*fpy);
	}
	public String toString(){
		return "Decomposition Method III";
	}
}