function plotOutput( out )

if strcmp(out.p.model(1:2),'FS')
    hf=figure(2);
    hf.Position(3)=1082;
    hf.Position(4)=420;
    
    subplot(1,3,1)
    plot(out.t,out.N)
    xlabel('t (s)')
    ylabel('N (cm^{-3})')

    subplot(1,3,2)
    plot(out.t,out.M_2)
    xlabel('t (s)')
    ylabel('M_2 (m^2 cm^{-3})')

    subplot(1,3,3)
    plot(out.t,out.M_3)
    xlabel('t (s)')
    ylabel('M_3 (m^3 cm^{-3})')

else
    hf=figure(2);
    hf.Position(3)=1082;
    hf.Position(4)=420;
    
    subplot(2,5,1)
    plot(out.t,out.Y(:,1))
    xlabel('t (s)')
    ylabel('N_{PL} (cm^{-3})')

    subplot(2,5,2)
    plot(out.t,out.alpha)
    xlabel('t (s)')
    ylabel('\alpha')
    ylim([out.alpha(end)-1 out.alpha(end)+1])

    subplot(2,5,3)
    plot(out.t,out.D2*1e9)
    xlabel('t (s)')
    ylabel('D_2 (nm)')
    
    subplot(2,5,4)
    plot(out.t,out.Y(:,2))
    xlabel('t (s)')
    ylabel('M_{PL,2} (m^2 cm^{-3})')
    
    subplot(2,5,5)
    plot(out.t,out.Y(:,3))
    xlabel('t (s)')
    ylabel('M_{PL,3} (m^3 cm^{-3})')

    subplot(2,5,6)
    plot(out.t,out.Y(:,4))
    xlabel('t (s)')
    ylabel('N_{LN} (cm^{-3})')

    subplot(2,5,7)
    plot(out.t,out.CMD*1e9)
    xlabel('t (s)')
    ylabel('CMD (nm)')

    subplot(2,5,8)
    plot(out.t,out.sigma)
    xlabel('t (s)')
    ylabel('\sigma')
    
    subplot(2,5,9)
    plot(out.t,out.Y(:,5))
    xlabel('t (s)')
    ylabel('M_{LN,2} (m^2 cm^{-3})')
    
    subplot(2,5,10)
    plot(out.t,out.Y(:,6))
    xlabel('t (s)')
    ylabel('M_{LN,3} (m^3 cm^{-3})')
end

end

