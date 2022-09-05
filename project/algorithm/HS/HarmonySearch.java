package project.algorithm.HS;

import project.data.*;
import project.files.SaveFile;
import project.start.*;

import java.util.*;

public class HarmonySearch 
{
	private int harmony_size;
	private int number_iterations;
	// Considering rate of harmony search, to determine the new harmony value
	private double harmony_search_rate;
	// Pitch adjusting rate: If the new harmony value is selected from the harmony memory then adjust to a neighbour value.
	private double pitch_adjusting_rate;
	private double neighbour_gap; 
	
	/* In screening problems, it is unlikely the optimisation will be able to obtain a set of cut offs that maximise both sensitivity and specificity. 
	 * It is highly likely that improvements in sensitivity will result in decreased specificity. 
	 * However when determining the pareto front of possible results, one can established the desired values for sensitivity and specificity to guide the search and to increase selection pressure towards 
	 * those individuals that obtain a pair (sensitivity, specificity) that is acceptable/desirable for that particular problem. This is why the variables below have been created, to steer the search.
	 */
	private double sensitivity_goal;
	private double specificity_goal;
	
	private ArrayList<Harmony> harmony_memory = null;
	private ArrayList<Harmony> best_harmony = null;
	private ArrayList<Harmony> update_harmony = null; 
	
	private int threshold;
	
	public HarmonySearch(Global_GA_HS_Configurations global_configurations) 
	{
		this.harmony_size = global_configurations.getHarmonySize();
		this.number_iterations = global_configurations.getHSIterations();
		this.harmony_search_rate = global_configurations.getHarmonySearchRate();
		this.pitch_adjusting_rate = global_configurations.getPitchAdjustingRate();
		this.neighbour_gap = global_configurations.getNeighbourGap();
		this.harmony_memory = new ArrayList<Harmony>();
		this.best_harmony = new ArrayList<Harmony>();
		this.update_harmony = new ArrayList<Harmony>();
		this.sensitivity_goal = global_configurations.getSensitivityGoal();
		this.specificity_goal = global_configurations.getSpecificityGoal();
		this.threshold = global_configurations.getThreshold();
	}
	
	public void run(int[] chromosome, ArrayList<Peptide> peptide_names, ArrayList<Sample> sample_values, Global_GA_HS_Configurations global_configurations, int iteration, int total_generations) 
	{

		initialHarmonyBuild(chromosome, peptide_names, sample_values, iteration, total_generations, global_configurations.random_generator);

		updateHarmony(chromosome, peptide_names, sample_values, global_configurations.random_generator, iteration, total_generations);
	}
	
	public void initialHarmonyBuild(int[] chromosome, ArrayList<Peptide> peptide_names, ArrayList<Sample> sample_values, int iteration, int total_iteration, RandomGenerator random_generator) 
	{
		for (int i = 0; i < this.harmony_size; i++) 
		{
			Harmony harmony = new Harmony(peptide_names.size()); // Harmony has the same size as the size of the chromosome in this implementation
			harmony.build(chromosome, peptide_names, random_generator);

			harmony.calculateVar(sample_values, this.threshold);
			harmony.calculateFitness(this.sensitivity_goal, this.specificity_goal, iteration, total_iteration);
			
			this.harmony_memory.add(harmony);
		}
		
		Collections.sort(this.harmony_memory, Harmony.sortByFITNESS);
		//PrintOnScreenBest();
	}
	
	public void updateHarmony(int[] chromosome, ArrayList<Peptide> peptide_names, ArrayList<Sample> sample_values, RandomGenerator random_generator, int iteration, int total_iteration) 
	{
		int iter = 0;

		while (iter < this.number_iterations) 
		{
			Harmony new_harmony = new Harmony(peptide_names.size());
			new_harmony.setChromosome(chromosome);

			for (int i = 0; i < chromosome.length; i++) 
			{
				if (chromosome[i] == 1) // this means the peptide has been included in the panel, and therefore a cut off needs to be established
				{
					double probability_harmony_search_rate = random_generator.getRandomValue(0, 1);
					double probability_pitch_adjusting_rate = random_generator.getRandomValue(0, 1);
					
					if (probability_harmony_search_rate > this.harmony_search_rate) // change the harmony to another feasible cut off within the window size
					{
						int index = random_generator.getRandomIndex(peptide_names.get(i).getCutOffsWindowSize());
						if(index != 0) // empty window. all peptides have same value
						{
							// exploration of the search space
							new_harmony.setHarmonyIndexValue(i, peptide_names.get(i).getCutOffValue(index));
						}
					} 
					else 
					{
						int index = (int)random_generator.getRandomIndex(this.harmony_size);
						double val = this.harmony_memory.get(index).getHarmonyIndexValue(i);
						
						if (probability_pitch_adjusting_rate < this.pitch_adjusting_rate && peptide_names.get(i).getCutOffsWindowSize() != 0 ) 
						{
							// This adjustment of the pitch means getting a neighbour cut off - exploitation of the search space locally
							Peptide peptide = peptide_names.get(i);
							int value_index = peptide.getCutOffValueIndex(val);
							
							if(value_index == -1) 
							{
								value_index = 0;
							}
							int neighbour_index = peptide.getNeighbourCutOffIndex(value_index, this.neighbour_gap);
							val = peptide.getCutOffValue(neighbour_index);
						}
						new_harmony.setHarmonyIndexValue(i, val);
					}
				} 
				else 
				{
					new_harmony.setHarmonyIndexValue(i, 0);
				}
			}
			
			// The existing hamony is replaced only if the offspring produced has a higher fitness value
			new_harmony.calculateVar(sample_values, threshold);
			new_harmony.calculateFitness(this.sensitivity_goal, this.specificity_goal, iteration, total_iteration);
				
			// version 1.2: minimising fitness
			if (new_harmony.getFitness() < this.harmony_memory.get(this.harmony_size-1).getFitness()) 
			{
				this.harmony_memory.remove(this.harmony_size-1);
				this.harmony_memory.add(new_harmony);
				
				if (new_harmony.getFitness() < this.harmony_memory.get(0).getFitness()) 
				{
					this.update_harmony.add(new_harmony);
				}
						
				Collections.sort(this.harmony_memory, Harmony.sortByFITNESS);
			}
			
			iter++;
			
			Harmony best_harmony = new Harmony(peptide_names.size());
			best_harmony = this.harmony_memory.get(0);
			this.best_harmony.add(best_harmony);
		}
	}
	
	public ArrayList<Harmony> getHarmonyMemory() 
	{
		return (ArrayList<Harmony>)this.harmony_memory.clone();
	}
	
	public Harmony getFinalBestHarmony() 
	{
		return this.harmony_memory.get(0);
	}
	
	public ArrayList<Harmony> getBestHarmony() 
	{
		return (ArrayList<Harmony>) this.best_harmony.clone();
	}
	
	public ArrayList<Harmony> getUpdateHarmony () 
	{
		return (ArrayList<Harmony>) this.update_harmony.clone(); 
	}
	
	public void PrintOnScreen() {
		for (int i = 0; i < this.harmony_memory.size(); i++) 
		{
			this.harmony_memory.get(i).PrintOnScreen();
		}
	}
	
	public void PrintOnScreenBest() 
	{
		this.harmony_memory.get(0).PrintOnScreen();
	}
}
