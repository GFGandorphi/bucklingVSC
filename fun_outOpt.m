function [stop,options,optchanged] = fun_outOpt(options,optimvalues,flag)
	
    stop = false;
    optchanged = false;
    if mod(optimvalues.iteration,10)==0 || flag=="done"
	    fprintf("  %s %4d : ",flag,optimvalues.iteration)
	    % fprintf(" [ <%5.1f|%5.1f> / <%5.1f|%5.1f> ]s ",optimvalues.x(1), optimvalues.x(2), optimvalues.x(3), optimvalues.x(4))
	    fprintf(" [ %5.1f %5.1f %5.1f %5.1f ]  ",optimvalues.x(1), optimvalues.x(2), optimvalues.x(2), optimvalues.x(1))
	    fprintf("  f = %5.1f ",optimvalues.fval)
	    % fprintf("  || best : [ <%5.1f|%5.1f> / <%5.1f|%5.1f> ]s ",optimvalues.bestx(1), optimvalues.bestx(2), optimvalues.bestx(3), optimvalues.bestx(4))
	    fprintf("  || best : [ %5.1f %5.1f %5.1f %5.1f ] ",optimvalues.bestx(1), optimvalues.bestx(2), optimvalues.bestx(2), optimvalues.bestx(1))
	    fprintf("  f = %5.1f  \n",optimvalues.bestfval)
    end
	
end