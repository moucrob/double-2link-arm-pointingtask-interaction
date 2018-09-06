function yesorno = japanpost(length,width,height)

if length > 150 %cm
    yesorno = 0 ; return ; end

if length + 2*(height + width) > 300
    yesorno = 0; %0 = not sendable
else
    yesorno = 1 ; end

end