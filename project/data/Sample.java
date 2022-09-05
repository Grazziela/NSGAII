package project.data;

import java.util.Arrays;

public class Sample 
{
	private int sample_id = 0;
	private String sample_name = "";
	private int number_peptides; // How many peptides/antigens are investigated in the database
	private boolean is_a_test_sample;
	private double[] Peptide_values = null;
	
	public Sample() {}
	
	public Sample(int number_peptides) 
	{
		this.number_peptides = number_peptides;
		this.Peptide_values = new double[number_peptides];
		this.is_a_test_sample = false;
	}
	
	public Sample(int sample_id, int number_peptides, boolean is_a_test_sample, double[] Peptide_values) {
		this.sample_id = sample_id;
		this.number_peptides = number_peptides;
		this.is_a_test_sample = is_a_test_sample;
		this.Peptide_values = new double[number_peptides];
		for (int i = 0; i < number_peptides; i++) {
			this.Peptide_values[i] = Peptide_values[i];
		}
	}
	
	public Sample(int sample_id, String[] input_data_file_values) 
	{
		this.sample_id = sample_id;
		this.sample_name = String.valueOf(sample_id);
		int header_length = input_data_file_values.length;
		this.number_peptides = header_length - 1;
		
		// In the input data file, the last column will be the class, which is a numerical value, i.e., test = 1 and control = 0	
		this.is_a_test_sample = (Integer.parseInt(input_data_file_values[header_length - 1]) == 0) ? false : true;

		this.Peptide_values = new double[this.number_peptides];
		
		for(int i = 0; i < this.number_peptides; i++) 
		{
			this.Peptide_values[i] = Double.parseDouble(input_data_file_values[i]);
		}
	}
	
	public String getSampleName(){
		return this.sample_name;
	}
	
	public int getPatientID() {
		return this.sample_id;
	}
	
	public int getNumberPeptides() {
		return this.number_peptides;
	}

	public boolean is_a_test_sample() {
		return this.is_a_test_sample;
	}

	public void setCancer(boolean is_a_test_sample) {
		this.is_a_test_sample = is_a_test_sample;
	}
	
	public double getPeptideValue(int index) {
		return this.Peptide_values[index];
	}
	
	public double[] getPeptideValues() {
		return this.Peptide_values.clone();
	}

	public void setPeptideValues(double[] Peptide_values) {
		int length = Peptide_values.length;
		this.Peptide_values = new double[length];
		for (int i = 0; i < length; i++) {
			this.Peptide_values[i] = Peptide_values[i];
		}
	}
	
	public void fixPeptideValue(int index, double value) {
		this.Peptide_values[index] = value;
	}
		
	public void printOnScreeen() {
		System.out.println(this.sample_id + ": "+ Arrays.toString(this.Peptide_values));
		System.out.println("The class of this sample is "+this.is_a_test_sample);
	}
	
	@Override
	public Sample clone() {
		Sample cloned = null;
		try {
			cloned = (Sample) super.clone();
		} 
		catch (CloneNotSupportedException e) 
		{			
		}		
		return cloned;
	}
}
