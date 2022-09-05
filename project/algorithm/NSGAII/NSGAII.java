package project.algorithm.NSGAII;

import project.start.RandomGenerator;
import project.data.*;
import project.algorithm.HS.*;
import project.files.*;

import java.util.*;

public class NSGAII {

	// NSGA-II parameters
	private int population_size;
	private int number_generations;
	private double crossover_rate;
	private double mutation_rate;
	private int chromosome_length;
	private int random_seed;
	private int panel_size;

	private ArrayList<Individual> initial_population = null;
	private ArrayList<Individual> current_population = null;
	private ArrayList<Individual> best_population = null;
	private ArrayList<Individual> updated_population = null;

	private double sensitivity_goal = 0;
	private double specificity_goal = 0;

	private int harmony_size_set = 0;

	private long start_time = 0;

	public NSGAII(Global_GA_HS_Configurations global_configurations, int chromossome_size) 
	{
		this.population_size = global_configurations.getPopulationSize();
		this.number_generations = global_configurations.getGAGenerations();
		this.crossover_rate = global_configurations.getCrossoverRate();
		this.mutation_rate = global_configurations.getMutationRate();
		this.chromosome_length = global_configurations.getChromosomeSize(); // The lenght of the chromossome is determined by the size of the peptides investigated (number of columns - 1)
		this.random_seed = global_configurations.getRandomSeed();
		this.panel_size = global_configurations.getPanelSize();

		this.initial_population = new ArrayList<Individual>();
		this.current_population = new ArrayList<Individual>();
		this.best_population = new ArrayList<Individual>();
		this.updated_population = new ArrayList<Individual>();

		this.sensitivity_goal = global_configurations.getSensitivityGoal();
		this.specificity_goal = global_configurations.getSpecificityGoal();
		this.harmony_size_set = global_configurations.getHarmonySize();
	}

	public void run(ArrayList<Peptide> peptide_names, ArrayList<Sample> sample_values, String save_folder, String final_result_folder, Global_GA_HS_Configurations global_configurations, long start_time) 
	{
		System.out.println("\n\nNSGA Started...\n");
		this.start_time = start_time;

		RandomGenerator rnd = new RandomGenerator(this.random_seed);
		System.out.println("\n\nChromosome Length: " + this.chromosome_length);

		PopulationInitialise(peptide_names, sample_values, save_folder, 0,  global_configurations);
		this.current_population = (ArrayList<Individual>)this.initial_population.clone();

		int front_size = 0;
		int generation = 0;
		while (generation < this.number_generations) 
		{
			System.out.println("\nGENERATION: " + generation);
			front_size = buildNextGeneration(chromosome_length, panel_size, peptide_names, sample_values, save_folder, generation, global_configurations);
			generation++;
		}

		SaveFile sf = new SaveFile();
		double running_time = getRunningTime();
		sf.saveNSGAIIResults(front_size, running_time, final_result_folder, "Pareto_Front", this.current_population, peptide_names, global_configurations);
	}

	/* In the population initialisation, the binary chromosome is created (1-hot indicators for presence or absence in the panel). 
	 * Subsequently, the optimal cut offs for those peptides are determined via harmony search.
	*/
	public void PopulationInitialise(ArrayList<Peptide> peptides, ArrayList<Sample> sample_values, String save_folder, int generation, Global_GA_HS_Configurations global_configurations) 
	{

		SaveFile sf = new SaveFile();


		for (int i = 0; i < this.population_size; i++) 
		{
			Individual chromosome = new Individual(chromosome_length);

			chromosome.buildIndividual(panel_size, global_configurations.random_generator, sample_values); //CMH Issue here, as the call to build the bit string resets the random seed call hence all bit strings on start are the same!
			System.out.println("\nBuilding Individual " + i);

			System.out.println("Peptides indexes: " + chromosome.getChomosomeIndex());
			
			HarmonySearch harmony_search = new HarmonySearch(global_configurations);

			// Optimising the cut offs for the peptides of the generated chromosome using Harmony Search
			harmony_search.run(chromosome.getChomosome(), peptides, sample_values, global_configurations, generation, this.number_generations);
			// This function also calculated the chromosome fitness
			chromosome.setHarmony(harmony_search.getFinalBestHarmony());
			
			this.initial_population.add(chromosome);
		}

		ArrayList<ArrayList<Integer>> fronts = new ArrayList<ArrayList<Integer>>();
		fronts = NonDominatedSort(this.initial_population);

		double running_time = getRunningTime();
		sf.saveNSGAIIResults(save_folder, "Initialisation", this.initial_population, running_time, peptides);
	}

	public int buildNextGeneration(int chromosome_length, int panel_size, ArrayList<Peptide> Peptides, ArrayList<Sample> Sample_info, String save_folder, int generation, Global_GA_HS_Configurations global_configurations)
	{
		RandomGenerator rnd = new RandomGenerator(this.random_seed);
		SaveFile sf = new SaveFile();

		ArrayList<Integer> population_index = new ArrayList<Integer>();
		for (int i = 0; i < this.population_size; i++) 
		{
			population_index.add(i);
		}


		ArrayList<Individual> overall_population = new ArrayList<Individual>();

		while (population_index.size() > 0) 
		{
			int parent1_index = selectParent(population_index, rnd);
			population_index.remove((Integer) parent1_index);
			int parent2_index = selectParent(population_index, rnd);
			population_index.remove((Integer) parent2_index);

			Individual parent1 = new Individual(chromosome_length);
			Individual parent2 = new Individual(chromosome_length);
			parent1 = this.current_population.get(parent1_index);
			parent2 = this.current_population.get(parent2_index);

			int[] parent1_chr = parent1.getChomosome();
			int[] parent2_chr = parent2.getChomosome();

			int[] child1_chr = new int[parent1_chr.length];
			int[] child2_chr = new int[parent2_chr.length];

			for (int i = 0; i < chromosome_length; i++) {
				if (rnd.getRandomValue(0, 1) > this.crossover_rate) 
				{
					child1_chr[i] = parent1_chr[i];
					child2_chr[i] = parent2_chr[i];
				} 
				else 
				{
					child1_chr[i] = parent2_chr[i];
					child2_chr[i] = parent1_chr[i];
				}

				if (rnd.getRandomValue(0, 1) < this.mutation_rate) 
				{
					child1_chr[i] = 1 - child1_chr[i];
				}

				if (rnd.getRandomValue(0, 1) < this.mutation_rate) 
				{
					child2_chr[i] = 1 - child2_chr[i];
				}
			}

			Individual child1 = new Individual(chromosome_length);
			buildIndividual(child1, child1_chr, Peptides, Sample_info, generation, global_configurations);
			overall_population.add(child1);

			Individual child2 = new Individual(chromosome_length);
			buildIndividual(child2, child2_chr, Peptides, Sample_info, generation, global_configurations);
			overall_population.add(child2);

			parent1.NSGAIIInitialisation();
			overall_population.add(parent1);
			parent2.NSGAIIInitialisation();
			overall_population.add(parent2);
		}

		ArrayList<ArrayList<Integer>> fronts = new ArrayList<ArrayList<Integer>>();
		fronts = NonDominatedSort(overall_population);

		ArrayList<Individual> aux_population  = new ArrayList<Individual>();
		int front_index = -1;
		
		// The population will be ranked. Fronts will indicate the individuals belonging to a certain rank, in order of fitness (rank 0 best fitness)
		// The new population will be the result of the best ranked individuals, based on their fitness and consequently rank.
		// For this reason, there is no need to worry about elitism being implemented here
		for (int i = 0; i < fronts.size(); i++) 
		{
			front_index = i;
			ArrayList<Integer> front = fronts.get(i);
			if (fronts.get(i).size() + aux_population.size() <= this.population_size) 
			{
				for (int j = 0; j < front.size(); j++) {
					Individual chr = new Individual(chromosome_length);
					chr = overall_population.get(front.get(j));
					aux_population.add(chr);
				}
			} 
			else 
			{
				break;
			}
		}

		// Gets the list of non-dominated solutions
		ArrayList<Integer> front = fronts.get(front_index);
		ArrayList<Individual> front_population = new ArrayList<Individual>();
		for (int i = 0; i < front.size(); i++) 
		{
			Individual chr = new Individual(chromosome_length);
			chr = overall_population.get(front.get(i));
			front_population.add(chr);
		}

		CrowdingDistanceAssignment(front_population);
		Collections.sort(front_population, Individual.sortByDISTANCE);
		
		for (int i = 0; i < front_population.size(); i++) 
		{
			if (aux_population.size() == this.population_size) 
			{
				break;
			} 
			else if (aux_population.size() < this.population_size) 
			{
				Individual chr = new Individual(chromosome_length);
				chr = front_population.get(i);
				aux_population.add(chr);
			}
		}

		this.current_population.clear();
		this.current_population = (ArrayList<Individual>)aux_population.clone();
		String iter = "generation" + generation;
		double running_time = getRunningTime();
		if (generation%10 == 0)
		{
			sf.saveNSGAIIResults(save_folder, iter, aux_population, running_time, Peptides);
		}

		return fronts.get(0).size();
	}

	public ArrayList<ArrayList<Integer>> NonDominatedSort(ArrayList<Individual> population) 
	{

		ArrayList<ArrayList<Integer>> fronts = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> first_front = new ArrayList<Integer>();
		ArrayList<Integer> uncheck_list = new ArrayList<Integer>();

		for (int i = 0; i < population.size(); i++) 
		{
			int number_dominated_by_individual = 0;
			ArrayList<Integer> dominators = new ArrayList<Integer>();
			//double sens = population.get(i).getSensitivity();
			//double spec = population.get(i).getSpecificity();
			double fitness = population.get(i).getFitnessValue();

			for (int j = 0; j < population.size(); j++) {
				if (j != i) 
				{
					//double sens_compare = population.get(j).getSensitivity();
					//double spec_compare = population.get(j).getSpecificity();
					double fitness_compare = population.get(j).getFitnessValue();

					// checking if there are better solutions. Original code based on sens/spec but I think this is wrong. Needs to be base on fitness
					//if ((sens <= sens_compare && spec <= spec_compare) && (sens < sens_compare || spec < spec_compare)) 
				    if (fitness > fitness_compare)
					{
						number_dominated_by_individual++;
					} 
					//else if((sens_compare <= sens && spec_compare <= spec) && (sens_compare < sens || spec_compare < spec)) 
				    else if(fitness < fitness_compare)
					{
						dominators.add(j);
					}
				}
			}

			population.get(i).setDominators(dominators);
			population.get(i).setNumberDominates(number_dominated_by_individual);
			uncheck_list.add(i);

			if (number_dominated_by_individual == 0) 
			{
				first_front.add(i);
				uncheck_list.remove((Integer) i);
				population.get(i).setRank(0);
			}
		}
		fronts.add(first_front);

		int rank = 1;
		while(!uncheck_list.isEmpty()) 
		{
			ArrayList<Integer> current_front = new ArrayList<Integer>();
			ArrayList<Integer> front = new ArrayList<Integer>();
			front = (ArrayList<Integer>)fronts.get(fronts.size()-1).clone();
			for (int i = 0; i < front.size(); i++) {
				ArrayList<Integer> domators = population.get(front.get(i)).getDominators();
				for (int j = 0; j < domators.size(); j++) {
					population.get(domators.get(j)).decreaseNumberDominates();
				}
			}

			for (int i = 0; i < uncheck_list.size(); i++) 
			{
				if (population.get(uncheck_list.get(i)).getNumberDominates() == 0) 
				{
					current_front.add(uncheck_list.get(i));
					population.get(uncheck_list.get(i)).setRank(rank);
					uncheck_list.remove((Integer) uncheck_list.get(i));
					i--;
				}
			}
			fronts.add(current_front);
			rank++;
		}

		return (ArrayList<ArrayList<Integer>>)fronts.clone();
	}

	public void CrowdingDistanceAssignment(ArrayList<Individual> chromosomes) 
	{
		double large_value = 10000.0;
		Collections.sort(chromosomes, Individual.sortBySENS);
		for (int i = 0; i < chromosomes.size(); i++) {
			double distance = chromosomes.get(i).getCrowdDistance();
			if (i == 0 || i == chromosomes.size()-1) {
				distance = distance + large_value;
			} else {
				distance = distance + (chromosomes.get(i-1).getSensitivity() - chromosomes.get(i+1).getSensitivity())/(chromosomes.get(0).getSensitivity() - chromosomes.get(chromosomes.size()-1).getSensitivity());
			}
			chromosomes.get(i).setCrowdDistance(distance);
		}

		Collections.sort(chromosomes, Individual.sortBySPEC);
		for (int i = 0; i < chromosomes.size(); i++) {
			double distance = chromosomes.get(i).getCrowdDistance();
			if (i == 0 || i == chromosomes.size()-1) {
				distance = distance + large_value;
			} else {
				distance = distance + (chromosomes.get(i-1).getSpecificity() - chromosomes.get(i+1).getSpecificity())/(chromosomes.get(0).getSpecificity() - chromosomes.get(chromosomes.size()-1).getSpecificity());
			}
			chromosomes.get(i).setCrowdDistance(distance);
		}
	}

	// Selection via tournament
	public int selectParent(ArrayList<Integer> population_index, RandomGenerator rnd) 
	{
		int selection = - 1;
		// Selecting for tournament the top 50% best individuals to increase the selection pressure
		int parent1_index = population_index.get(rnd.getRandomIndex(population_index.size()/2));
		int parent2_index = population_index.get(rnd.getRandomIndex(population_index.size())/2);

		// Tournament selection
		if (this.current_population.get(parent1_index).getRank() <= this.current_population.get(parent2_index).getRank()) 
		{
			selection = parent1_index;
		} 
		else 
		{
			selection = parent2_index;
		}

		return selection;
	}

	public void buildIndividual(Individual chr, int[] chromosome_index, ArrayList<Peptide> Peptide_info, ArrayList<Sample> Sample_info, int generation, Global_GA_HS_Configurations global_configurations)
	{
		chr.buildIndividual(chromosome_index);
		if (chr.getTotalGen() > panel_size)
		{
			chr.flipGen(panel_size, this.random_seed);
		} 
		else if(chr.getTotalGen() < panel_size)
		{
			chr.antiflipGen(panel_size, this.random_seed);
		}
		HarmonySearch hs = new HarmonySearch(global_configurations);
		hs.run(chr.getChomosome(), Peptide_info, Sample_info, global_configurations, generation, this.number_generations);
		chr.setHarmony(hs.getFinalBestHarmony()); // Fitness is calculated here
	}


	public double getRunningTime() 
	{
		long end_time = System.currentTimeMillis();
		double running_time = (end_time - this.start_time)*1.0/1000; 
		return running_time; 
	}


	public int getPopulationSize() 
	{
		return this.population_size;
	}

	public int getNumberGeneration() {
		return this.number_generations;
	}

	public double getCrossoverRate() {
		return this.crossover_rate;
	}

	public double getMutationRate() {
		return this.mutation_rate;
	}

	public int getChromosomeLength() {
		return this.chromosome_length;
	}

	public ArrayList<Individual> getInitialPopulation() {
		return (ArrayList<Individual>) this.initial_population.clone();
	}

	public ArrayList<Individual> getBestPopulation() {
		return (ArrayList<Individual>) this.best_population.clone();
	}

	public ArrayList<Individual> getUpdatedPopulation() {
		return (ArrayList<Individual>) this.updated_population.clone();
	}

	// After the population is ordered (ascending), the best individuals will be in the position zero
	public double[] getBestValues() {
		double[] values = new double[2];
		values[0] = this.best_population.get(0).getSensitivity();
		values[1] = this.best_population.get(0).getSpecificity();
		return values;
	}
}
