function [ im_enhanced,enhancement_lut,best_fitness,pheromone_map,best_chromosome,fitness_per_iteration,last_enhancing_part,elapsed_time ] = imenhance( im_inputt,iteration_num,sa_disable )
%[Output Image,TF,Best Fitness,Pheromone Map,Best Chromosome,Fitness per Iteration,
%Last Enhancing Part,Elapsed Time]=imenhance(Input Image, Iteration Numbers, 'no_sa')
%Function takes input image, number of interations & a string which
%controls the SA part to be disabled or not. And gives enhanced image,
%transfer function's look up table, best fitness obtained, final pheromone
%map, best chromosome that finally assigned to all ants and fitness per 
%iteration diagram's look up table. Also an string (Last Enhancing Part) 
%specifies, between ACO or SA, which of them lastly edited the output 
%transfer function. Final output is total time elapsed during process, too.
%Image Enhances by Using ACO, GA and SA. Ant Colony Optimization Used to
%Find Best Transfer Function of Intensities. Genetic Algorithm Used for
%Optimize ACO. Simulated Annealing Used for Local Search After ACO's Global
%Search.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Preparing Part

%start timing
tic;

%global declarations
global im_input min_input_intensity max_input_intensity min_input_intensity_plus_1 max_input_intensity_minus_1 relative_input_range;
im_input=im_inputt;

%set random stream
get_clock=clock;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(25000000*get_clock(4:6))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Input-Output Arguments Checking Part

%set default value of iterationnum
if nargin==1
    iteration_num=100;
    sa_disable='sa';    %keeps SA enabled
elseif nargin==2
    sa_disable='sa';    %keeps SA enabled
end

%round iterationnum
iteration_num=round(iteration_num);

%check number of inputs & outputs
narginchk(1,3);
nargoutchk(1,8);

%check correction of inputs class
if ~isa(im_input,'uint8')
    error('Input image must be gray scale');
end

if ~isa(iteration_num,'double')
    error('Number of iterations must be double');
end

if ~ischar(sa_disable)
    error('Disabling string of SA must be char');
end

%check that SA must be enabled or not
sa_en=~strcmp(sa_disable,'no_sa');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Input Image Processing

%extract min & max intensity of original image
image_range=double(im2uint8(stretchlim(im_input)))+1;
min_input_intensity=image_range(1);
max_input_intensity=image_range(2);
min_input_intensity_plus_1=min_input_intensity+1;
max_input_intensity_minus_1=max_input_intensity-1;
relative_input_range=(max_input_intensity-min_input_intensity)/255;

%show original image
%figure,imshow(im_input);
%title('Original Image');
%drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Fitness Function Estimation Part

    function [fitness,im_current]=fitnesscalc(lut_current)
        
        im_current=intlut(im_input,lut_current);
        %<<<<<<REMOVED>>>>>> co_occurrence=graycoprops(graycomatrix(im_current),'contrast');
        fitness=(std2(im_current)*entropy(im_current)*(mean2(abs(imfilter(im_current,fspecial('sobel')))+abs(imfilter(im_current,fspecial('sobel')')))))^(1/3); %<<<<<<REMOVED>>>>>> (mean2(stdfilt(im_current))),(co_occurrence.Contrast);
        
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Main Body

%database for interrupt main body
sa_schedule=[10:10:round(40*iteration_num/100),round(46*iteration_num/100):6:round(70*iteration_num/100),round(73*iteration_num/100):3:iteration_num];
sa_point_num=[2*ones(1,ceil(4*iteration_num/100)),4*ones(1,ceil(5*iteration_num/100)),6*ones(1,ceil(10*iteration_num/100))];
sa_ant_num=[1*ones(1,ceil(4*iteration_num/100)),2*ones(1,ceil(5*iteration_num/100)),4*ones(1,ceil(10*iteration_num/100))];
sa_duration=[3*ones(1,ceil(4*iteration_num/100)),6*ones(1,ceil(5*iteration_num/100)),12*ones(1,ceil(10*iteration_num/100))];
ga_schedule=[5:5:round(25*iteration_num/100),round(31*iteration_num/100):6:round(55*iteration_num/100),round(62*iteration_num/100):7:round(90*iteration_num/100)];
sa_interrupt_num=1;
ga_interrupt_num=1;

%set initial value of some variables
pheromone_map=zeros(256);
best_fitness=0;
gene_fitness=zeros(1,10);
genetic_elite=zeros(1,10);
iterations_untill_ga=0;
best_ga_fitness=0;
enhancement_lut=zeros(1,256);
elite_pheromone_trace=zeros(256);
fitness_per_iteration=zeros(1,iteration_num);
last_enhancing_part='Ant Colony Optimization';

%set initial temperature
temperature=200;

%initialize ga
probability_factor_vector_up(2:2:20)=3*rand(1,10);
probability_factor_vector_right(2:2:20)=3*rand(1,10);
impact_factor_vector_alpha(2:2:20)=5*rand(1,10);
impact_factor_vector_beta(2:2:20)=5*rand(1,10);
routing_factor_vector(2:2:20)=randi(150,1,10);
probability_factor_vector_up(1:2:20)=probability_factor_vector_up(2:2:20);
probability_factor_vector_right(1:2:20)=probability_factor_vector_right(2:2:20);
impact_factor_vector_alpha(1:2:20)=impact_factor_vector_alpha(2:2:20);
impact_factor_vector_beta(1:2:20)=impact_factor_vector_beta(2:2:20);
routing_factor_vector(1:2:20)=routing_factor_vector(2:2:20);

%main loop
for iteration=1:iteration_num
    
    %call ant colony function
    [fitness_vector,lut_matrix,pheromone_trace_matrix,enhancement_lut,best_fitness,elite_pheromone_trace,last_enhancing_part]=aco(elite_pheromone_trace,pheromone_map,enhancement_lut,best_fitness,probability_factor_vector_up,probability_factor_vector_right,impact_factor_vector_alpha,impact_factor_vector_beta,100+routing_factor_vector,last_enhancing_part);
    
    %sum of fitness for each gene
    iterations_untill_ga=iterations_untill_ga+1;
    gene_fitness=gene_fitness+[fitness_vector(1)+fitness_vector(2) fitness_vector(3)+fitness_vector(4) fitness_vector(5)+fitness_vector(6) fitness_vector(7)+fitness_vector(8) fitness_vector(9)+fitness_vector(10) fitness_vector(11)+fitness_vector(12) fitness_vector(13)+fitness_vector(14) fitness_vector(15)+fitness_vector(16) fitness_vector(17)+fitness_vector(18) fitness_vector(19)+fitness_vector(20)];
    genetic_elite=max([max(fitness_vector(1),fitness_vector(2)) max(fitness_vector(3),fitness_vector(4)) max(fitness_vector(5),fitness_vector(6)) max(fitness_vector(7),fitness_vector(8)) max(fitness_vector(9),fitness_vector(10)) max(fitness_vector(11),fitness_vector(12)) max(fitness_vector(13),fitness_vector(14)) max(fitness_vector(15),fitness_vector(16)) max(fitness_vector(17),fitness_vector(18)) max(fitness_vector(19),fitness_vector(20))],genetic_elite);
    
    %check for sa turn
    if iteration==sa_schedule(sa_interrupt_num) && sa_en
        
        %call simulated annealing function
        [enhancement_lut,best_fitness,pheromone_trace_matrix,elite_pheromone_trace,fitness_vector,last_enhancing_part]=sa(temperature,enhancement_lut,best_fitness,pheromone_trace_matrix,elite_pheromone_trace,fitness_vector,lut_matrix,sa_point_num(sa_interrupt_num),sa_ant_num(sa_interrupt_num),sa_duration(sa_interrupt_num),last_enhancing_part);
        
        %preparing next turn check for sa
        sa_interrupt_num=sa_interrupt_num+1;
        
        %reduce temperature
        temperature=temperature*(1-0.5*300/iteration_num);
        
        %prevent running sa in next iterations if it is no planning to run
        %it
        if sa_interrupt_num>numel(sa_schedule)
            sa_interrupt_num=sa_interrupt_num-1;
        end
        
    end
    
    %check for ga turn
    if iteration==ga_schedule(ga_interrupt_num)
        
        %calculate fitness of chromosomes
        ga_fitness=(gene_fitness/iterations_untill_ga)+genetic_elite;
        
        %save best chromosome and its fitness
        [best_ga_fitness_candidate,best_ga_index_candidate]=max(ga_fitness);
        if best_ga_fitness_candidate>=best_ga_fitness
            best_ga_fitness=best_ga_fitness_candidate;
            best_chromosome=[probability_factor_vector_up(2*best_ga_index_candidate),probability_factor_vector_right(2*best_ga_index_candidate),impact_factor_vector_alpha(2*best_ga_index_candidate),impact_factor_vector_beta(2*best_ga_index_candidate),routing_factor_vector(2*best_ga_index_candidate)];
        end
        
        %call genetic algorithm function
        [probability_factor_vector_up,probability_factor_vector_right,impact_factor_vector_alpha,impact_factor_vector_beta,routing_factor_vector]=ga(probability_factor_vector_up,probability_factor_vector_right,impact_factor_vector_alpha,impact_factor_vector_beta,routing_factor_vector,ga_fitness);
        
        %preparing next turn check for ga
        ga_interrupt_num=ga_interrupt_num+1;
        
        %check for last ga turn
        if ga_interrupt_num>numel(ga_schedule)
            ga_interrupt_num=ga_interrupt_num-1; %prevents running ga in next iterations
            %assign best gene to all ants
            probability_factor_vector_up(:)=best_chromosome(1);
            probability_factor_vector_right(:)=best_chromosome(2);
            impact_factor_vector_alpha(:)=best_chromosome(3);
            impact_factor_vector_beta(:)=best_chromosome(4);
            routing_factor_vector(:)=best_chromosome(5);
        end
        
        %reset fitness summers
        iterations_untill_ga=0;
        gene_fitness=zeros(1,10);
        genetic_elite=zeros(1,10);

    end
    
    %write pheromone map (fitness/(30*best fitness) adding)
    %evaporation_rate=0.4
    %plus elitism section
    weighted_pheromone_trace_matrix=(repmat(permute(fitness_vector,[3 1 2]),[256 256 1])).*pheromone_trace_matrix;
    pheromone_map=0.6*pheromone_map+(sum(weighted_pheromone_trace_matrix,3)/(30*best_fitness))+(0.1*elite_pheromone_trace/30); %0.1 is quarter evaporation rate
           
    %printing results in process time
    fitness_per_iteration(iteration)=best_fitness;
  %  if mod(iteration,50)==0
  %      [best_fitness,im_enhanced]=fitnesscalc(enhancement_lut);
  %      figure,imshow(im_enhanced);
  %      title({'Enhanced Image by ',last_enhancing_part});
  %      figure,plot(enhancement_lut,'color',[0.3 0.1 0.1],'linewidth',3);
  %      axis tight;
  %      title('Transfer Function');
  %      figure,imshow(pheromone_map,[]);
  %      title('Pheromone Map');
  %      figure,plot(fitness_per_iteration(1:iteration));
  %      title('Fitness per Iteration');
  %      drawnow;
  %      disp('******');
  %      disp('Iteration=');
  %      disp(iteration);
  %      disp('Best Chromosome=');
  %      disp(best_chromosome+[0 0 0 0 100]);
  %      toc;
  %      disp('******');
  %  end

end

%create output image
[best_fitness,im_enhanced]=fitnesscalc(enhancement_lut);

%create other outputs
best_chromosome(5)=best_chromosome(5)+100;
elapsed_time=toc;

%compare with fitness of original image
input_lut=uint8(0:1:255);
input_fitness=fitnesscalc(input_lut);
if best_fitness<input_fitness
    replace=input('Best found fitness is lower than fitness of original image. To show original image as result press "o" or press any other key to show processed image: ','s');
    if strcmp(replace,'o') || strcmp(replace,'O')
        im_enhanced=im_input;
        enhancement_lut=input_lut;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Genetic Algorithm (GA) Part

    function [probability_factor_vector_up,probability_factor_vector_right,impact_factor_vector_alpha,impact_factor_vector_beta,routing_factor_vector]=ga(probability_factor_vector_up,probability_factor_vector_right,impact_factor_vector_alpha,impact_factor_vector_beta,routing_factor_vector,ga_fitness)
        
        %parent selection (roulette wheel)
        [sorted_ga_fitness,ga_fitness_index]=sort(ga_fitness);
        selected_chromosome_1=ga_fitness_index((cumsum(sorted_ga_fitness)/sum(sorted_ga_fitness))>=rand());
        censored_ga_fitness=ga_fitness;
        censored_ga_fitness(selected_chromosome_1(1))=0;    %prevents twice selection of a chromosome
        [sorted_censored_ga_fitness,censored_ga_fitness_index]=sort(censored_ga_fitness);
        selected_chromosome_2=censored_ga_fitness_index((cumsum(sorted_censored_ga_fitness)/sum(sorted_censored_ga_fitness))>=rand());
        
        %sort parents
        parent_index=[selected_chromosome_1(1) selected_chromosome_2(1)];
        [~,sorted_parent_index]=sort(ga_fitness(parent_index));
        
        %worst chromosome selection
        worst_chromosome_index=ga_fitness_index(1);
        
        %change worst chromosome if it is weaker parent too
        if parent_index(sorted_parent_index(1))==ga_fitness_index(1)
            worst_chromosome_index=ga_fitness_index(2);
        end
                
        %uniform crossover
        child_1=[probability_factor_vector_up(2*parent_index(1)),probability_factor_vector_right(2*parent_index(1)),impact_factor_vector_alpha(2*parent_index(1)),impact_factor_vector_beta(2*parent_index(1)),routing_factor_vector(2*parent_index(1))];
        child_2=[probability_factor_vector_up(2*parent_index(2)),probability_factor_vector_right(2*parent_index(2)),impact_factor_vector_alpha(2*parent_index(2)),impact_factor_vector_beta(2*parent_index(2)),routing_factor_vector(2*parent_index(2))];
        if rand<=0.85
            for gene_counter=1:5
                if rand<0.5
                    crossover_stash=child_1(gene_counter);
                    child_1(gene_counter)=child_2(gene_counter);
                    child_2(gene_counter)=crossover_stash;
                end
            end
        end
        
        %one point mutation
        %offspring 1
        if rand<0.05
            mutation_gene=randi(5);
            if mutation_gene==5
                child_1(mutation_gene)=child_1(mutation_gene)+randi(30)-15;
                if child_1(mutation_gene)>150 %keep genes in the permitted range
                    child_1(mutation_gene)=150;
                end
            elseif mutation_gene>=3
                child_1(mutation_gene)=child_1(mutation_gene)+rand-0.5;
                if child_1(mutation_gene)>5 %keep genes in the permitted range
                    child_1(mutation_gene)=5;
                end
            elseif mutation_gene<=2
                child_1(mutation_gene)=child_1(mutation_gene)+(0.4*rand)-0.2;
                if child_1(mutation_gene)>2 %keep genes in the permitted range
                    child_1(mutation_gene)=2;
                end
            end
            %keep genes in the permitted range
            if child_1(mutation_gene)<0
                child_1(mutation_gene)=0;
            end
        end
        %offspring 2
        if rand<0.05
            mutation_gene=randi(5);
            if mutation_gene==5
                child_2(mutation_gene)=child_2(mutation_gene)+randi(30)-15;
                if child_2(mutation_gene)>150 %keep genes in the permitted range
                    child_2(mutation_gene)=150;
                end
            elseif mutation_gene>=3
                child_2(mutation_gene)=child_2(mutation_gene)+rand-0.5;
                if child_2(mutation_gene)>5 %keep genes in the permitted range
                    child_2(mutation_gene)=5;
                end
            elseif mutation_gene<=2
                child_2(mutation_gene)=child_2(mutation_gene)+(0.4*rand)-0.2;
                if child_2(mutation_gene)>2 %keep genes in the permitted range
                    child_2(mutation_gene)=2;
                end
            end
            %keep genes in the permitted range
            if child_2(mutation_gene)<0
                child_2(mutation_gene)=0;
            end
        end
        
        %replacing weaker parent by offspring 1
        probability_factor_vector_up(2*parent_index(sorted_parent_index(1))-1:2*parent_index(sorted_parent_index(1)))=child_1(1);
        probability_factor_vector_right(2*parent_index(sorted_parent_index(1))-1:2*parent_index(sorted_parent_index(1)))=child_1(2);
        impact_factor_vector_alpha(2*parent_index(sorted_parent_index(1))-1:2*parent_index(sorted_parent_index(1)))=child_1(3);
        impact_factor_vector_beta(2*parent_index(sorted_parent_index(1))-1:2*parent_index(sorted_parent_index(1)))=child_1(4);
        routing_factor_vector(2*parent_index(sorted_parent_index(1))-1:2*parent_index(sorted_parent_index(1)))=child_1(5);
        
        %replacing worst chromosome by offspring 2
        probability_factor_vector_up(2*worst_chromosome_index-1:2*worst_chromosome_index)=child_2(1);
        probability_factor_vector_right(2*worst_chromosome_index-1:2*worst_chromosome_index)=child_2(2);
        impact_factor_vector_alpha(2*worst_chromosome_index-1:2*worst_chromosome_index)=child_2(3);
        impact_factor_vector_beta(2*worst_chromosome_index-1:2*worst_chromosome_index)=child_2(4);
        routing_factor_vector(2*worst_chromosome_index-1:2*worst_chromosome_index)=child_2(5);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Ant Colony Optimization (ACO) Part

    function [fitness_vector,lut_matrix,pheromone_trace_matrix,enhancement_lut,best_fitness,elite_pheromone_trace,last_enhancing_part]=aco(elite_pheromone_trace,pheromone_map_in,enhancement_lut,best_fitness,probability_factor_vector_up,probability_factor_vector_right,impact_factor_vector_alpha,impact_factor_vector_beta,routing_factor_vector,last_enhancing_part)
        
        %lut and pheromone map generating through walking ants
        lut_matrix=zeros(20,256,'uint8');
        pheromone_trace_matrix=zeros(256,256,20);
        fitness_vector=zeros(1,20);
        routing_factor_vector_up=relative_input_range*routing_factor_vector;
        
        for ant=1:20
            intensity=1;
            height=1;
            current_lut=zeros(1,256,'uint8');
            pheromone_trace=zeros(256);
            pheromone_trace(256,1)=1;   %each ant has a trace in (1,1)
            
           while ~(intensity==256 && height==256)
               
               %boundary checking
               if intensity<min_input_intensity || intensity>max_input_intensity || height==256
                   next_point_probability=[0 0 1];
               elseif intensity==max_input_intensity
                   next_point_probability=[1 0 0];
               else
                   %next point probability calculation ([up upper-right right])
                   next_point_probability=([((pheromone_map_in(257-(height+1),intensity)+1)^impact_factor_vector_alpha(ant))*((probability_factor_vector_up(ant)*(1+((intensity-min_input_intensity)/routing_factor_vector_up(ant))^10))^(impact_factor_vector_beta(ant))),((pheromone_map_in(257-(height+1),intensity+1)+1)^impact_factor_vector_alpha(ant)),((pheromone_map_in(257-(height),intensity+1)+1)^impact_factor_vector_alpha(ant))*((probability_factor_vector_right(ant)*(1+(height/routing_factor_vector(ant))^10))^(impact_factor_vector_beta(ant)))])/(((pheromone_map_in(257-(height+1),intensity)+1)^impact_factor_vector_alpha(ant))*((probability_factor_vector_up(ant)*(1+((intensity-min_input_intensity)/routing_factor_vector_up(ant))^10))^(impact_factor_vector_beta(ant)))+((pheromone_map_in(257-(height+1),intensity+1)+1)^impact_factor_vector_alpha(ant))+((pheromone_map_in(257-(height),intensity+1)+1)^impact_factor_vector_alpha(ant))*((probability_factor_vector_right(ant)*(1+(height/routing_factor_vector(ant))^10))^(impact_factor_vector_beta(ant))));
               end
               
               %selection (roulette wheel)
               [sorted_probability,probability_index]=sort(next_point_probability);
               selected_directions=probability_index(cumsum(sorted_probability)>=rand());
               
               %moving ant
               switch selected_directions(1)
                   case 1   %up
                       height=height+1;
                   case 2   %upper-right
                       height=height+1;
                       intensity=intensity+1;
                   case 3   %right
                       intensity=intensity+1;
               end
               
               %write into lut
               current_lut(intensity)=height-1;
               
               %record ant's pheromone trace
               pheromone_trace(257-height,intensity)=1;
               
           end
           
           %calculation of fitness
           fitness_vector(ant)=fitnesscalc(current_lut);
        
           %finding best lut
           if fitness_vector(ant)>=best_fitness
               enhancement_lut=current_lut;
               best_fitness=fitness_vector(ant);
               elite_pheromone_trace=pheromone_trace;   %save elite pheromone trace
               last_enhancing_part='Ant Colony Optimization';
           end
           
           %save current pheromone trace
           pheromone_trace_matrix(:,:,ant)=pheromone_trace;
           
           %save current lut
           lut_matrix(ant,:)=current_lut;
           
        end
                        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Simulated Annealing (SA) Part

    function [enhancement_lut,best_fitness,pheromone_trace_matrix,elite_pheromone_trace,fitness_vector,last_enhancing_part]=sa(temperature,enhancement_lut,best_fitness,pheromone_trace_matrix,elite_pheromone_trace,fitness_vector,lut_matrix,sa_point_num,sa_ant_num,sa_duration,last_enhancing_part)
        
        %prepare for optimizing best lut
        lut_matrix(21,:)=enhancement_lut;
        selected_lut=21;
        fitness_vector(21)=best_fitness;
        pheromone_trace_matrix(:,:,21)=elite_pheromone_trace;
        
        lower_bound=0;
        
        for lut_num=1:sa_ant_num+1
                        
            last_point=min_input_intensity-1;
            
            %obtaining the items to be optimized
            optimization_lut=lut_matrix(selected_lut,:);
            optimization_fitness=fitness_vector(selected_lut);
            optimization_pheromone_trace=pheromone_trace_matrix(:,:,selected_lut);
            
            for point_num=1:sa_point_num
                
                %selecting a point to be optimized
                last_point=last_point+randi(max_input_intensity-sa_point_num+point_num-last_point);
                active_point=last_point;
                                
                for sa_cycle=1:sa_duration
                    
                    %prepare for next iteration
                    new_lut=optimization_lut;
                    new_pheromone_trace=optimization_pheromone_trace;
                    new_active_point=active_point;
                    break_search=1;
                    
                    %reset probability of neighbors
                    neighborhood_probability=[1,1,1,1,0,1,1,1,1];
                    
                    if active_point==min_input_intensity
                        
                        %restriction of neighborhood at vertical border of lut
                        neighborhood_probability(1:3:7)=0;
                        
                        %restriction of neighborhood to keep lut
                        %monotonically increasing
                        if optimization_lut(active_point)==optimization_lut(active_point+1)
                            neighborhood_probability(2)=0;
                        end
                        
                        if optimization_lut(active_point)==optimization_lut(active_point+2)
                            neighborhood_probability(3)=0;
                        end
                        
                    elseif active_point==max_input_intensity
                        
                        %restriction of neighborhood at vertical border of lut
                        neighborhood_probability(3:3:9)=0;
                        
                        %restriction of neighborhood to keep lut
                        %monotonically increasing
                        if optimization_lut(active_point)==optimization_lut(active_point-1)
                            neighborhood_probability(8)=0;
                        end
                        
                        if optimization_lut(active_point)==optimization_lut(active_point-2)
                            neighborhood_probability(7)=0;
                        end
                    
                    elseif active_point==min_input_intensity_plus_1
                        
                        %restriction of neighborhood to keep lut
                        %monotonically increasing
                        if optimization_lut(active_point)==optimization_lut(active_point+1)
                            neighborhood_probability(1:2)=0;
                        end
                                        
                        if optimization_lut(active_point)==optimization_lut(active_point-1)
                            neighborhood_probability(8:9)=0;
                        end
                    
                        if optimization_lut(active_point)==optimization_lut(active_point+2)
                            neighborhood_probability(3)=0;
                        end
                        
                    elseif active_point==max_input_intensity_minus_1
                        
                        %restriction of neighborhood to keep lut monotonically increasing
                        if optimization_lut(active_point)==optimization_lut(active_point+1)
                            neighborhood_probability(1:2)=0;
                        end
                        
                        if optimization_lut(active_point)==optimization_lut(active_point-1)
                            neighborhood_probability(8:9)=0;
                        end
                        
                        if optimization_lut(active_point)==optimization_lut(active_point-2)
                            neighborhood_probability(7)=0;
                        end
                        
                    else
                        
                        %restriction of neighborhood to keep lut monotonically increasing
                        if optimization_lut(active_point)==optimization_lut(active_point+1)
                            neighborhood_probability(1:2)=0;
                        end
                                        
                        if optimization_lut(active_point)==optimization_lut(active_point-1)
                            neighborhood_probability(8:9)=0;
                        end
                    
                        if optimization_lut(active_point)==optimization_lut(active_point+2)
                            neighborhood_probability(3)=0;
                        end
                    
                        if optimization_lut(active_point)==optimization_lut(active_point-2)
                            neighborhood_probability(7)=0;
                        end
                        
                    end
                    
                    %restriction of neighborhood at horizontal border of lut
                    if optimization_lut(active_point)==0
                        neighborhood_probability(7:9)=0;
                    elseif optimization_lut(active_point)==255
                        neighborhood_probability(1:3)=0;
                    end
                    
                    %move to new point
                    neighbor_num=numel(find(neighborhood_probability));
                    for neighbor_counter=1:neighbor_num
                        
                        %randomly suggest a direction to move
                        remaining_neighbor=find(neighborhood_probability);
                        suggested_direction=remaining_neighbor(randi(numel(remaining_neighbor)));
                        
                        %remove suggested direction from remaining
                        %directions
                        neighborhood_probability(suggested_direction)=0;
                        
                        %calculate effect of suggested direction
                        if suggested_direction==1
                            new_lut(active_point-1)=optimization_lut(active_point)+1;
                            new_lut(active_point)=optimization_lut(active_point)+1;
                            new_pheromone_trace(1:256,active_point)=0;  %clean ant's previous way
                            new_pheromone_trace(double(256-new_lut(active_point)),active_point)=1;  %write ant's new way
                            if active_point~=max_input_intensity
                                new_pheromone_trace(double(257-new_lut(active_point)),active_point+1)=0;  %clean ant's previous way
                            end
                            new_pheromone_trace(double(256-new_lut(active_point)):double(256-optimization_lut(active_point-1)),active_point-1)=1;  %write ant's new way                            
                        elseif suggested_direction==2
                            new_lut(active_point)=optimization_lut(active_point)+1;
                            new_pheromone_trace(double(256-new_lut(active_point)),active_point)=1;  %write ant's new way
                            if active_point~=max_input_intensity
                                new_pheromone_trace(double(257-new_lut(active_point)),active_point+1)=0;  %clean ant's previous way
                            end
                        elseif suggested_direction==3
                            new_lut(active_point+1)=optimization_lut(active_point)+1;
                            new_active_point=active_point+1;
                            new_pheromone_trace(double(256-optimization_lut(new_active_point)):double(255-new_lut(new_active_point)),new_active_point)=0;   %clean ant's previous way
                            if new_active_point~=max_input_intensity
                                new_pheromone_trace(double(256-optimization_lut(new_active_point+1)):double(256-new_lut(new_active_point)),new_active_point+1)=1;  %write ant's new way
                            end
                            new_pheromone_trace(double(256-new_lut(new_active_point)),new_active_point)=1;  %write ant's new way
                        elseif suggested_direction==4
                            new_lut(active_point-1)=optimization_lut(active_point);
                            new_active_point=active_point-1;
                            new_pheromone_trace(double(256-new_lut(new_active_point)):256-optimization_lut(new_active_point),new_active_point)=1;  %write ant's new way
                            new_pheromone_trace(double(257-new_lut(new_active_point)):256,active_point)=0;  %clean ant's previous way
                        elseif suggested_direction==6
                            new_lut(active_point+1)=optimization_lut(active_point);
                            new_active_point=active_point+1;
                            new_pheromone_trace(1:256,new_active_point)=0;   %clean ant's previous way
                            if new_active_point~=max_input_intensity
                                new_pheromone_trace(double(256-optimization_lut(new_active_point+1)):double(256-new_lut(new_active_point)),new_active_point+1)=1;  %write ant's new way
                            end
                            new_pheromone_trace(double(256-new_lut(new_active_point)),new_active_point)=1;  %write ant's new way
                        elseif suggested_direction==7
                            new_lut(active_point-1)=optimization_lut(active_point)-1;
                            new_active_point=active_point-1;
                            new_pheromone_trace(double(256-new_lut(new_active_point)):256-optimization_lut(new_active_point),new_active_point)=1;  %write ant's new way
                            new_pheromone_trace(double(257-new_lut(new_active_point)):256,active_point)=0;  %clean ant's previous way
                            new_pheromone_trace(double(255-new_lut(new_active_point)),new_active_point)=0;  %clean ant's previous way
                        elseif suggested_direction==8
                            new_lut(active_point)=optimization_lut(active_point)-1;
                            new_pheromone_trace(double(256-new_lut(active_point)),active_point)=1;  %write ant's new way
                            if active_point~=max_input_intensity
                                new_pheromone_trace(double(255-new_lut(active_point)),active_point+1)=1;  %write ant's new way
                            end
                            new_pheromone_trace(double(255-new_lut(active_point)),active_point)=0;  %clean ant's previous way
                        elseif suggested_direction==9
                            new_lut(active_point+1)=optimization_lut(active_point)-1;
                            new_lut(active_point)=optimization_lut(active_point)-1;
                            new_pheromone_trace(double(256-new_lut(active_point)),active_point)=1;  %write ant's new way
                            new_pheromone_trace(double(255-new_lut(active_point)),active_point)=0;  %clean ant's previous way
                            new_pheromone_trace(1:256,active_point+1)=0;   %clean ant's previous way
                            if active_point~=max_input_intensity_minus_1
                                new_pheromone_trace(double(256-optimization_lut(active_point+2)):double(256-new_lut(active_point)),active_point+2)=1;  %write ant's new way
                            end
                            new_pheromone_trace(double(256-new_lut(active_point)),active_point+1)=1;  %write ant's new way
                        end
                        
                        %decide to select the suggested direction or not
                        new_fitness=fitnesscalc(new_lut);
                        if new_fitness>=optimization_fitness || rand<=exp((new_fitness-optimization_fitness)/(0.05*optimization_fitness)*temperature)
                            optimization_lut=new_lut;
                            optimization_pheromone_trace=new_pheromone_trace;
                            active_point=new_active_point;
                            optimization_fitness=new_fitness;
                            break_search=0;
                            break
                        else
                            new_active_point=active_point;
                            new_lut=optimization_lut;
                            new_pheromone_trace=optimization_pheromone_trace;
                        end
                        
                    end
                    
                    %check for break search or not
                    if break_search==1
                        break
                    end
                    
                    %search for best fitness
                    if optimization_fitness>=best_fitness
                        best_fitness=optimization_fitness;
                        enhancement_lut=optimization_lut;
                        elite_pheromone_trace=optimization_pheromone_trace;
                        last_enhancing_part='Simulated Annealing';
                    end
                                    
                end
                
            end
            
            %saving the optimized items
            fitness_vector(selected_lut)=optimization_fitness;
            pheromone_trace_matrix(:,:,selected_lut)=optimization_pheromone_trace;
            
            %selecting one lut to be optimized
            selected_lut=lower_bound+randi(20/sa_ant_num);
            lower_bound=lower_bound+(20/sa_ant_num);
            
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end