
class PolyElement implements Func{
	protected Double coefficient;
	protected int exponent;
	
	PolyElement(Double co, int ex){
		this.coefficient = co;
		this.exponent = ex;
	}
	
	PolyElement(String input){
		//expect format like "3:6", or "C:E"
		String[] temp = input.split(":", 2);
		this.coefficient = Double.parseDouble(temp[0]);
		this.exponent = Integer.parseInt(temp[1]);
	}
	
	public String toString(){
		if(exponent == 0){return coefficient.toString();}
		if(exponent == 1){return coefficient +"*x";}
		if(coefficient == 0){return "";}
		
		return coefficient +"*x^" + exponent;
	}
	
	public String toString(String type){
		if(!type.equals("latex")){
			return this.toString();
		}
		
		if(exponent == 0){return coefficient.toString();}
		if(exponent == 1){return coefficient +"\\cdot x";}
		if(coefficient == 0){return "";}
		
		return coefficient +"\\cdot x^{" + exponent+"}";
	}
	
	public Double evaluate(Double x){
		return this.coefficient * Math.pow(x, this.exponent);
	}
	
	public Func differentiate(){
		return new PolyElement(this.coefficient * this.exponent ,this.exponent-1);
	}
	
}