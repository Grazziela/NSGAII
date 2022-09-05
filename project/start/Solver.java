package project.start;

import project.algorithm.*;
import project.algorithm.HS.*;
import project.algorithm.NSGAII.*;
import project.data.*;
import project.files.*;

import java.util.*;

public class Solver 
{
	private ArrayList<Sample> sample_values = null;
	private ArrayList<Peptide> peptide_names = null;
	
	private int panel_size;

	private double[] best_values = null;
	
	private long start_time = 0;
	
	public Solver(int panel_size) 
	{
		this.sample_values = new ArrayList<Sample>();
		this.peptide_names = new ArrayList<Peptide>();
		this.panel_size = panel_size;
		
	}
	
	// Method that creates the GA and populations
	public void start(ReadInputDatasetFile data, String save_folder, String final_result_folder, Global_GA_HS_Configurations global_configurations) 
	{

		this.start_time = System.currentTimeMillis();
		this.best_values = new double[2];

		System.out.println("\nThe start running time: " + this.start_time + "    Panel size:" + this.panel_size + "    Threshold: " + global_configurations.getThreshold());
		
		this.sample_values = data.getSampleValues();	
		this.peptide_names = data.getPeptideNames();
				
		System.out.println("Sensitivity Goal: " + global_configurations.getSensitivityGoal() + "    Specificity Goal: " + global_configurations.getSpecificityGoal());
		
		NSGAII nsga = new NSGAII(global_configurations, this.peptide_names.size());
		
		System.out.println("\n\nPopulation Size: " + nsga.getPopulationSize());		
		nsga.run(this.peptide_names, this.sample_values, save_folder, final_result_folder, global_configurations, this.start_time);		
	}

	
	public int[] getIndexOfPeptide(String[] panel) {
		int[] chromosome_index = new int[this.panel_size];
		for (int i = 0; i < panel.length; i++) {
			for (int j = 0; j < this.peptide_names.size(); j++) {
				if (this.peptide_names.get(j).getPeptideName().equals(panel[i])) {
					chromosome_index[i] = j;
					break;
				}
			}
		}
		
		return chromosome_index.clone();
	}
	
	public double[] getBestValues() {
		return this.best_values;
	}
}
