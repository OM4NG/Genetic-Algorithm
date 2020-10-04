%Code for Genetic Algorithm using Roulette Wheel
%Made by: Omang Saxena 194103424
%Part of Soft computing Course (Assignment-2)
%{
N = population size
max_val = max. value that the variable can take
min_val = min. value that the variable can take
%L = length of binary string for each variable
max_gen = max. number of generations generated
pc = probability of cross-over
pm = probability of mutation
population = random generated population
x = population at each iteration
decode_x = decimal value of binary string
power = to kee count of opwer while converting binary to decimal
value_for_fun = value of variables to be used in the function
fitness_val = fitness_val calculated using value_for_fun
temporary = temp. variable for performing two point cross-over
o = stores original funcion value. (to be minimised)
p1 and p2 =  parents
ch1 and ch2 = children
point1 and point2 = points for cross-over
cum_prob = to store cumulative probability
%}
%%
clc;
N = input("Enter population size: ");
min_val = input("Enter minimum value: ");
max_val = input("Enter maximum value: ");
L = input("Enter length of each element in binary string: ");
max_gen = input("Enter maximum number of generations: ");
pc = input("Enter cross-over probability: ");
pm = input("Enter mutation probability: ");
%%
gen = 1;

%Generating random population
population = randi([0 1],N,2*L);

%To store new population
x = zeros(N,2*L,max_gen);

for i = 1:N
    x(i,1:L,1) = population(i,1:L);
    x(i, L+1:2*L,1) = population(i,L+1:2*L);
end

%To store decoded values
decode_x = zeros(N,2);

%% Starting the GA
while(gen<max_gen)
    for i = 1:N
        %calculating the decoded values
        num1 =0;
        num2 =0;
    
        power = L-1;
        for k = 1:L
            num1 = num1 + x(i,k,gen)*2^power;
            power = power -1;
        end
            decode_x(i,1) = num1; %decoded value 1
   
        power = L-1;
        for k = L+1:2*L
            num2  = num2  + x(i,k,gen)*2^power;
            power = power -1;
        end
        decode_x(i,2) = num2; %decoded value 2
    end

    %% Calculating Values to be substituted in the function for fitness calculation
    val_for_fun = zeros(N,2);

    for i = 1:N
        val_for_fun(i,1) = min_val + (max_val - min_val)/(2^L - 1)*decode_x(i,1);
        val_for_fun(i,2) = min_val + (max_val - min_val)/(2^L - 1)*decode_x(i,2);
    end

    %% Storing the different fitness values
    fitness_val = zeros(N);
    for i = 1:N
        fitness_val(i) = func(val_for_fun(i,1),val_for_fun(i,2));
    end

    fitness_prob = fitness_val/sum(fitness_val);
    %% Calculating the cumulative probability for Roulette wheel selection
    cum_prob = zeros(N);
    for i =1:N
        if (i == 1)
            cum_prob(i) = fitness_prob(i);
        else
            for k = 1:i
                cum_prob(i) = cum_prob(i) + fitness_prob(k);
            end
        end
    end
    
    
    %% Defining the mating pool
    mating_pool = zeros(N,2*L);

    %Randomly copying entries in the mating pool
    for i = 1:N
        p = randi([0.0 1.0],1,1);
        if (p<cum_prob(1))
            p = 1;
        elseif (p>cum_prob(1) && p<cum_prob(2))
            p = 2;
        elseif(p>cum_prob(2) && p<cum_prob(3))
            p = 3;
        elseif(p>cum_prob(3) && p<cum_prob(4))
            p = 4;
        else
            p = 5;
        end
        mating_pool(i,:) = x(p,:,gen);
    end

    %% Making mating pairs
    mating_pair = zeros(floor(N/2), 2);
    if (mod(gen,2) ==0)
        for i = 1:floor(N/2)
            mating_pair(i,1) = i;
            mating_pair(i,2) = N + 1-i;
        end
    else
        for i = 1:floor(N/2)
            mating_pair(i,1) = i;
            mating_pair(i,2) = 2*i;
        end
    end
    
%% Performin the crossover between the mating pairs obtained
    for i = 1:floor(N/2)
        pp1 = mating_pair(i,1);
        pp2 = mating_pair(i,2);
        p1 = x(pp1,:,gen);
        p2 = x(pp2,:,gen);
        ch1 = p1;
        ch2 = p2;
    
        point1 = randi([1 2*L-1],1,1);
        point2 = randi([1 2*L-1],1,1);
        
        %This loop will help in getting two different values of point for
        %crossover
        while(point1==point2)
            point1 = randi([1 L-1],1,1);
            point2 = randi([1 L-1],1,1);
        end
    
        start = min(point1, point2) + 1;
        stop = max(point1,point2);
    
        %% Performing two-point cross-over
        if(randi([0.0 1.0],1,1)<pc)
            temporary  = ch1(1,start:stop);
            ch1(1,start:stop) = ch2(1,start:stop);
            ch2(1,start:stop) = temporary;
        end
    
        %% Performing mutation bit by bit
        for bit = 1:L
            if(randi([0.0 1.0],1,1)<pm)
                if(ch1(bit) ==0)
                    ch1(bit) = 1;
                else 
                    ch1(bit) = 0;
                end
            
                if(ch2(bit) ==0)
                    ch2(bit) = 1;
                else 
                    ch2(bit) = 0;
                end
            end
        end
        
        %% Storing the values in next generation 
    x(2*i -1,:,gen+1) = ch1;
    x(2*i,:,gen+1) = ch2;
    end
    gen = gen+1;
end
%% GA finished

%% Calculating the final values obtained 
    for i = 1:N
        %calculating the decoded values
        num1 =0;
        num2 =0;
    
        power = L-1;
        for k = 1:L
            num1 = num1 + x(i,k,max_gen)*2^power;
            power = power -1;
        end
            decode_x(i,1) = num1; %decoded value 1
   
        power = L-1;
        for k = L+1:2*L
            num2  = num2  + x(i,k,max_gen)*2^power;
            power = power -1;
        end
        decode_x(i,2) = num2; %decoded value 2
    end

        val_for_fun = zeros(N,2);
%% Calculating the value for function
    for i = 1:N
        val_for_fun(i,1) = min_val + (max_val - min_val)/(2^L - 1)*decode_x(i,1);
        val_for_fun(i,2) = min_val + (max_val - min_val)/(2^L - 1)*decode_x(i,2);
    end
    
    x(:,:,1);
    x(:,:,max_gen);
    fprintf("--------------------------------------------------\n");
    fprintf("  x1          x2           Max_func      Orig_func\n");
    
    %Storing the different fitness values
    fitness_val = zeros(N);
    
    for i = 1:N
        fitness_val(i) = func(val_for_fun(i,1),val_for_fun(i,2));
        
        % Variable 'o' stores the original function value which is to be
        % minimised
        o = original(val_for_fun(i,1),val_for_fun(i,2));
                
        fprintf("%.3f           %.3f            %.3f        %.3f\n",val_for_fun(i,1),val_for_fun(i,2),fitness_val(i),o);
    end
    

    



     
    
    
        
    
    