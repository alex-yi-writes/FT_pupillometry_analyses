function index1 = f_discontinuties(d,c)

    stop = d;
    start_1 = 0;
    index = [];
    
    list1 = stop(1:end-1);
    list2 = stop(2:end);
    s = 1; %start
    e = 1; %end

    for i=1:length(stop)-1
        if list2(i)-1 == list1(i) 
            e = e + 1;            
        elseif list2(i)-1 ~= list1(i) 
            index(end+1,:) = [stop(s),stop(e)];
            s = e+1;
            e = s;
        end
    
    end
    if i==1
       start_1 = i;
       index(end+1,:) = [stop(s),stop(e)];
    else
       index(end+1,:) = [stop(s),stop(e)];
    end


%% merge if two ranges are very close, say c=15 ms
    gap = index;
    e1 = 0;
    s1 = 0;
    j = 1;
    index1 = [];

    for i=1:size(gap,1)-1
        if e1 == 0 &&s1 == 0
           s1 = gap(i,1);
           e1 = gap(i,2);
        end
    
        if gap(i+1,1) - gap(i,2) <= c
             e1 = gap(i+1,2);
        end
        if gap(i+1,1) - gap(i,2)> c
              index1(j,1) = s1;
              index1(j,2) = e1;
              j=j+1;
              if i < size(gap,1)
                  s1 = gap(i+1,1);
                  e1 = gap(i+1,2);                  
              end
        end
        if i == length(gap)-1
           index1(j,1) = s1;
           index1(j,2) = e1;
        end
    end

end