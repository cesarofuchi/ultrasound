%%
function fdata=zerosFilter(data)

    n=size(data,2)
    dentro_intervalo=0;
    for i=1:n
        % valor dentro do intervalo?
        if(dentro_intervalo==0)
             if(data(i)==0)
                 if(i>1)                
                    amp1=data(i-1);
                    if(amp1>0)
                        dentro_intervalo=1;
                        i1=(i-1);
                        zcont=0;
                    end
                 end

            end

        else
            zcont=zcont+1;
            if(data(i)>0)
                dentro_intervalo=0;
                amp2=data(i);
                i2=i;
                %interpolar os dados
                data(i1:i2)=interpolar_dados([amp1 amp2], zcont);

            end

        end

    end
    fdata=data;
    %figure
    %plot(fdata)
    %grid on
end

function novo_vetor=interpolar_dados(vetor,zeros)
    x_desejado=1:1/(zeros+1):2;
    novo_vetor=interp1(1:2,vetor,x_desejado);
end
