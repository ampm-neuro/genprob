function ms=flood_fill(matrix,yc,xc)
% flood fill, scan line algoritm
% ms=flood_fill_new(I,r,tol)
% ms - pixels numbers that flooded
% numbering: number=y+matrix_size_Y*(x-1), x y - pixels coordinaties

% (yc, xc)- first point of selection

[matrix_size_Y, matrix_size_X]=size(matrix); % image size

% stak, where seeds will be stored
maximum_stack_size=10000; %needs to be higher for crazier shapes
stak = zeros(maximum_stack_size, 2, 'int32');
stak(1,1) = xc;
stak(1,2) = yc; %first seed
remaining_stack = 1; %starting stack length

% establish pixel store
maximum_field_size = 1000000; % margin
field_pixles=zeros(maximum_field_size,1,'int32'); % predifined array to increase speed
field_size=0; % 0 points initially

% to pixel number %??? %number=y+matrix_size_Y*(x-1), x y - pixels coordinaties
tn=@(xx,yy) yy+matrix_size_Y*(xx-1);


while true
    % get seed from stack:
    x_seed = stak(remaining_stack,1); %x coord of seed
    y_seed = stak(remaining_stack,2); %y coord of seed
    remaining_stack = remaining_stack-1;
    
    
    % LINE SCAN TO RIGHT
    
    %These variables alternate below true and false below, and their
    %current state is used to gate code behavior
    sku=false; % seed key, true if seed added, up
    skd=false; % same for down
    sku1=false; % this key is need to prewnt extra seed when move to left, up
    skd1=false; % same for left
    
    %xtt - active pixle within the row
    for xtt = x_seed:matrix_size_X %from x coord of seed to end of line
        if matrix(y_seed,xtt) == 1 %if current pixle needs to be changed
            % add pixels to field
            field_size=field_size+1;
            field_pixles(field_size)=tn(xtt,y_seed);
        else
            break;
        end
        
        % TRY TO ADD (just one) SEED UP
        if y_seed~=matrix_size_Y %if we're not already at the top
            if matrix(y_seed+1,xtt) == 1
                if ~sku %if we don't already have a seed from this line

                    if all(tn(xtt,y_seed+1)~=field_pixles(1:field_size)) % if free space (none of the above pixles are already in the field?)
                        % add to stack
                        remaining_stack=remaining_stack+1;
                        stak(remaining_stack,1)=xtt;
                        stak(remaining_stack,2)=y_seed+1;
                        sku=true; %now we have a seed
                    end
                end
            else
                sku=false; %don't have a seed
            end
            if xtt==x_seed
                sku1=sku; % memorize (that we're not already at the top?), will be used when to left
            end
        end
        
        % TRY TO ADD (just one) SEED DOWN
        if y_seed~=1 %if we're not already at the bottom
            if matrix(y_seed-1,xtt) == 1
                if ~skd
                    if all(tn(xtt,y_seed-1)~=field_pixles(1:field_size)) % if free space
                        % add to stack
                        remaining_stack=remaining_stack+1;
                        stak(remaining_stack,1)=xtt;
                        stak(remaining_stack,2)=y_seed-1;
                        skd=true;
                    end
                end
            else
                skd=false;
            end
            if xtt==x_seed
                skd1=skd; % memorize, will be used when to left
            end
        end
    end
    
    
    
    % LINE SCAN TO LEFT
    %sku=false; % seed key, true if seed added
    %skd=false;
    sku=sku1;
    skd=skd1;
    if x_seed~=1
        for xtt=(x_seed-1):-1:1 

            if matrix(y_seed,xtt) == 1 %if current pixle needs to be changed
                % add pixel
                field_size=field_size+1;
                field_pixles(field_size)=tn(xtt,y_seed);
            else
                break;
            end

            % TRY TO ADD (just one) SEED UP
            if y_seed~=matrix_size_Y
                if matrix(y_seed+1,xtt) == 1
                    if ~sku
                        if all(tn(xtt,y_seed+1)~=field_pixles(1:field_size)) % if free space
                            % add to stack
                            remaining_stack=remaining_stack+1;
                            stak(remaining_stack,1)=xtt;
                            stak(remaining_stack,2)=y_seed+1;
                            sku=true;
                        end
                    end
                else
                    sku=false;
                end
            end

            % TRY TO ADD (just one) SEED DOWN
            if y_seed~=1
                if matrix(y_seed-1,xtt) == 1
                    if ~skd
                        if all(tn(xtt,y_seed-1)~=field_pixles(1:field_size)) % if free space
                            % add to stack
                            remaining_stack=remaining_stack+1;
                            stak(remaining_stack,1)=xtt;
                            stak(remaining_stack,2)=y_seed-1;
                            skd=true;
                        end
                    end
                else
                    skd=false;
                end
            end
        end
    end
    
    if remaining_stack==0 % no more seed
        break; % stop
    end
    
    
end

ms=field_pixles(1:field_size);