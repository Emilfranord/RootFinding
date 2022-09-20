
class Poly implements Func{
	protected PolyElement[] elements;
	
	public Poly(String input){
		input = input.replace("^x", ":");
		input = input.replace("x", ":1");
		String[] splitInp = input.split(" ");
		PolyElement[] temp = new PolyElement[splitInp.length];
		
		int index = 0;
		for (String q : splitInp){
			temp[index] = new PolyElement(q);
			index++;
		}
		this.elements = temp;
		
	}
	
	public Poly(PolyElement[] input){
		this.elements = input;
	}
	
	public Double evaluate(Double x){
		Double sum = 0.0;
		for(PolyElement q : this.elements){
			sum += q.evaluate(x);
		}
		return sum;
	}
	
	public Func differentiate(){
		PolyElement[] future = new PolyElement[elements.length];
		for(int i= 0;i <elements.length;i++){
			future[i] = (PolyElement) elements[i].differentiate();
		} 
		return new Poly(future);
	}
	
	public String toString(){
		StringBuilder temp = new StringBuilder();
		boolean first = true;
		for(PolyElement q :elements){
			if (q.coefficient >=0){
				if (!first){
					temp.append("+");
				}
			}
			temp.append(q.toString());
			first = false;
		}
		return temp.toString();
		
	}
	
	public String toString(String type){
		if(!type.equals("latex")){
			return this.toString();
		}
		StringBuilder temp = new StringBuilder();
		boolean first = true;
		for(PolyElement q :elements){
			if (q.coefficient >=0){
				if (!first){
					temp.append("+");
				}
			}
			temp.append(q.toString("latex"));
			first = false;
		}
		return temp.toString();
	}
}