package project.data;

import java.util.*;

public class Peptide 
{

	private String Peptide_name;
	private int Peptide_id;
	private double max_value;
	private double min_value;
	private ArrayList<Double> peptide_values= null;

	private ArrayList<Double> cutoff_windows = null;

	public Peptide(String Peptide_name, int Peptide_id) 
	{
		this.Peptide_name = Peptide_name;
		this.Peptide_id = Peptide_id;
		initialise();
	}

	/*
	 * When establishing cut offs for Peptides, valid cut offs need to be located within the interval of values present in the dataset.
	 * Initialising the cut offs within those intervals avoids the creation of invalid chromossomes
	 */
	public void initialise() 
	{
		this.max_value = Double.MIN_VALUE;
		this.min_value = Double.MAX_VALUE;
		this.peptide_values = new ArrayList<Double>();
		this.cutoff_windows = new ArrayList<Double>();	
	}

	// This function makes sure the values attributed to the cut offs are valid and within the intervals present in the dataset
	public void update(double value) 
	{
		if (this.min_value > value) 
		{
			this.min_value = value;
		}

		if(this.max_value < value) 
		{
			this.max_value = value;
		}

		this.peptide_values.add(value);
	}

	/* To make the search and optimisation of cutt offs more efficient, we get the intervals between data points and generate a finite list of possible
	 * cut offs for that particular peptide/antigen. This method generates this list, which will be used to determine the values in the harmony of 
	 * the harmony search
	 */
	public void PreDefiningFeasibleCutOffs() 
	{
		Collections.sort(this.peptide_values);

		for (int i = 0; i < this.peptide_values.size() - 1; i++) 
		{
			double recent = this.peptide_values.get(i);
			double next = this.peptide_values.get(i+1);

			if (recent != next) 
			{
				// The cut off between two subsequent points is their midway value
				double window = recent + (next-recent)/2;
				this.cutoff_windows.add(window);
			}
		}
	}

	public int getNeighbourCutOffIndex(int index, double neighbor_range) 
	{
		int neighbour_index = -1;
		int change = (int) Math.round(this.cutoff_windows.size()*neighbor_range);

		if (index + change > this.cutoff_windows.size() && index - change > -1) 
		{
			neighbour_index = index - change;
		} 
		else if (index + change < this.cutoff_windows.size()) 
		{
			neighbour_index = index + change;
		} 
		else 
		{
			neighbour_index = index - change;
		}
		//System.out.println("neighbour index " + neighbour_index);
		if (neighbour_index < 0)
		{
			neighbour_index = 0;
		}

		return neighbour_index;
	}

	public int getCutOffValueIndex(double value) {
		int index = -1;
		int neighbor_index = -1;
		
		for (int i = 0; i < this.cutoff_windows.size(); i++) 
		{
			if (this.cutoff_windows.get(i) == value) 
			{
				index = i;
				break;
			} 
			else if(this.cutoff_windows.get(i) > value) 
			{
				neighbor_index = i-1;
				break;
			}
		}

		if (index == -1) {
			if (neighbor_index == -1) {
				index = this.cutoff_windows.size()-1;
			} else {
				index = neighbor_index;
			}
		}

		return index;
	}

	public double getNeighbourValue(double value) 
	{
		int neighbour_index = 0;
		return this.cutoff_windows.get(neighbour_index);
	}

	public String getPeptideName() 
	{
		return this.Peptide_name;
	}

	public void setPeptideName(String Peptide_name) 
	{
		this.Peptide_name = Peptide_name;
	}

	public int getPeptideId() 
	{
		return this.Peptide_id;
	}

	public void setPeptideId(int Peptide_id) 
	{
		this.Peptide_id = Peptide_id;
	}

	public double getMaxValue() 
	{
		return this.max_value;
	}

	public void setMaxValue(double max_value) {
		this.max_value = max_value;
	}

	public double getMinValue() {
		return this.min_value;
	}

	public void setMinValue(double min_value) {
		this.min_value = min_value;
	}

	public void setPeptideValues(ArrayList<Double> Peptide_values) 
	{
		Collections.copy(this.peptide_values, Peptide_values);
	}

	public ArrayList<Double> getPeptideValues() 
	{
		return this.peptide_values;
	}

	public ArrayList<Double> getWholeCutoffWindows() 
	{
		return this.cutoff_windows;
	}

	public int getCutOffsWindowSize() 
	{
		return this.cutoff_windows.size();
	}

	public double getCutOffValue(int index) 
	{
		return this.cutoff_windows.get(index);
	}

	@Override
	public Peptide clone() 
	{
		Peptide cloned = null;
		try {
			cloned = (Peptide) super.clone();
		} catch (CloneNotSupportedException e) {

		}

		return cloned;
	}
}
