function [] = pdc_xplot_title(alpha,metric)

alphastr = int2str(100*alpha);

   switch metric
      case 'euc'
         suptitle(['PDC (' '{\alpha = ' alphastr '%}' ')'])
      case 'diag'
         suptitle(['gPDC (' '{\alpha = ' alphastr '%}' ')'])
      case 'info'
         suptitle(['iPDC (' '{\alpha = ' alphastr '%}' ')'])
      otherwise
         error('Unknown metric.')
   end;