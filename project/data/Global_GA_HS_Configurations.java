package project.data;

import project.start.RandomGenerator;

/* 
 * This class works as a repository of global variables, which are used as parameters for the GA, HS and overall parameters of the solver.
 * The advantage is to have all parameters for a run centralised, which makes it easier to fine tune the solver for a particular problem.
 */
public class Global_GA_HS_Configurations {
	
	// Problem General Configurations
    private int panel_size; // size of the panel that will perform the screening (separation between tests and controls)
    
	/* This variable indicates how many Peptides/antigens etc need to have values above the cut off(s) to classify a sample as positive.
	* Default value for screening is 1. So for this case, if for one peptide within the panel the cut off established is below the value measured for the sample,
	* the sample classification is 1 (which means that the sample is from the test/cancer positive/etc. category)
	*/
    private int threshold; 
    
    // Target Configurations - these are the predefined values for sensitivity and specificity that will steer the search for panels within the Pareto Front
    private double sensitivity_goal;
    private double specificity_goal;
    
    private int random_seed; // To be used for random number generation for GA and HS
    private boolean generate_paretofront;
    
    // Genetic Algorithm parameters
    private int chromosome_size;
    private int population_size;
    private int GA_generations;
    private double crossover_rate;
    private double mutation_rate;
    
    // Harmony Search parameters
    private int harmony_size;
    private int harmony_search_iterations;
    private double harmony_search_rate;
    private double pitch_adjusting_rate;
    private double neighbour_gap;
    
    // Input data file parameters
    private int experiment_id;
    private String experiment_description;
    
    public RandomGenerator random_generator;
    
	public Global_GA_HS_Configurations(int panel_size, int chromosome_size, String experiment_description) 
	{			
		this.panel_size = panel_size;
		this.threshold = 1;
		
		this.sensitivity_goal = 0.60;
		this.specificity_goal = 1.0;
		this.random_seed = 123;
		
		this.chromosome_size = chromosome_size;
		this.population_size = 400;
		this.GA_generations = 180;
		this.crossover_rate = 0.65;
		this.mutation_rate = 0.15;
		
		this.harmony_size = 50;
		this.harmony_search_iterations = 200;
		this.harmony_search_rate = 0.6;
		this.pitch_adjusting_rate = 0.25;
		this.neighbour_gap = 0.05;
		
		this.experiment_description = experiment_description;
		
		this.generate_paretofront = true;
		
		 this.random_generator = new RandomGenerator(random_seed);
	}
	
	// Set/Get methods for the class attributes
    
    public int getPanelSize() 
    {
        return this.panel_size;
    }

    public void setPanelSize(int panel_size) 
    {
        this.panel_size = panel_size;
    }

    public int getThreshold() 
    {
        return this.threshold;
    }

    public void setThreshold(int threshold) 
    {
        this.threshold = threshold;
    }

    public double getSensitivityGoal() 
    {
        return this.sensitivity_goal;
    }

    public void setSensitivityGoal(double sensitivity_goal) 
    {
        this.sensitivity_goal = sensitivity_goal;
    }

    public double getSpecificityGoal() 
    {
        return this.specificity_goal;
    }

    public int getRandomSeed() 
    
    {
        return this.random_seed;
    }
    
    
    public boolean isGenerateParetoFront() 
    {
        return this.generate_paretofront;
    }

    public int getChromosomeSize() 
    {
        return this.chromosome_size;
    }

    
    public int getPopulationSize() 
    {
        return this.population_size;
    }


    public int getGAGenerations() 
    {
        return this.GA_generations;
    }


    public double getCrossoverRate() 
    {
        return this.crossover_rate;
    }

     public double getMutationRate() 
    {
        return this.mutation_rate;
    }

    public int getHarmonySize() 
    {
        return this.harmony_size;
    }

    public int getHSIterations() 
    {
        return this.harmony_search_iterations;
    }
 
    public double getHarmonySearchRate()
    {
        return this.harmony_search_rate;
    }
 
    public double getPitchAdjustingRate() 
    {
        return this.pitch_adjusting_rate;
    }

    public double getNeighbourGap() 
    {
        return this.neighbour_gap;
    }
    
     public String getExperimentDescription() 
    {
    	return this.experiment_description;
    }
}
