
1;   % penanda bahwa file ini adalah script (bukan file function)

% ------------------------- FUNGSI BANTU ----------------------------------
function y = cos_taylor(z, n)
  % pendekatan cos(z) dengan deret Taylor sampai suku ke-n
  y = zeros(size(z));
  for k = 0:n
    y = y + (-1)^k * z.^(2*k) / factorial(2*k);
  end
end

function y = exp_taylor(z, n)
  % pendekatan e^z dengan deret Taylor sampai suku ke-n
  y = zeros(size(z));
  for k = 0:n
    y = y + z.^k / factorial(k);
  end
end

% ------------------------- PROGRAM UTAMA ---------------------------------
clc; clear; close all;

% ---- data soal ----
lambda = 12;
v      = 36;
x      = 10;
h0     = 5;
N      = [23 24 25];            % batas iterasi
t      = linspace(0, 5, 501);   % 0 <= t <= 5
tol    = 1;                     % toleransi error = 1 %

% ---- 1) perhitungan eksak ----
h_eksak = h0 * (cos(t*v/lambda) + exp(x));

% ---- 2) perhitungan pendekatan (deret Taylor) ----
h_aprox = zeros(numel(N), numel(t));
for m = 1:numel(N)
  h_aprox(m,:) = h0 * (cos_taylor(t*v/lambda, N(m)) + exp_taylor(x, N(m)));
end

% ---- 3) galat ----
rerata_mutlak  = zeros(1, numel(N));
rerata_relatif = zeros(1, numel(N));
maks_relatif   = zeros(1, numel(N));
for m = 1:numel(N)
  galat_mutlak  = abs(h_eksak - h_aprox(m,:));
  galat_relatif = galat_mutlak ./ abs(h_eksak) * 100;     % dalam persen
  rerata_mutlak(m)  = mean(galat_mutlak);
  rerata_relatif(m) = mean(galat_relatif);
  maks_relatif(m)   = max(galat_relatif);
end

% ---- 4) tampilkan hasil ----
fprintf('==================================================================\n');
fprintf(' KASUS 1 : GALAT (ERROR PEMOTONGAN)\n');
fprintf('==================================================================\n');
fprintf(' lambda = %g, v = %g, x = %g, h0 = %g, 0 <= t <= 5 (%d titik)\n\n', ...
        lambda, v, x, h0, numel(t));

fprintf(' Nilai e^x pada x = %g\n', x);
fprintf('   eksak          : %.6f\n', exp(x));
for m = 1:numel(N)
  fprintf('   deret n = %-2d   : %.6f   (selisih %.6f)\n', ...
          N(m), exp_taylor(x, N(m)), exp(x) - exp_taylor(x, N(m)));
end

fprintf('\n h(t) pada t = 0 dan t = 5\n');
fprintf('   %-10s %18s %18s\n', 'Metode', 'h(t=0)', 'h(t=5)');
fprintf('   %-10s %18.6f %18.6f\n', 'Eksak', h_eksak(1), h_eksak(end));
for m = 1:numel(N)
  fprintf('   %-10s %18.6f %18.6f\n', sprintf('n = %d', N(m)), ...
          h_aprox(m,1), h_aprox(m,end));
end

fprintf('\n Perbandingan galat terhadap hasil eksak\n');
fprintf('   %-8s %20s %20s %20s\n', 'n', 'Rerata galat mutlak', ...
        'Rerata galat relatif', 'Galat relatif maks');
for m = 1:numel(N)
  fprintf('   %-8d %20.6f %17.6f %% %17.6f %%\n', N(m), ...
          rerata_mutlak(m), rerata_relatif(m), maks_relatif(m));
end

fprintf('\n Toleransi error = %g %%\n', tol);
for m = 1:numel(N)
  if rerata_relatif(m) < tol
    status = 'LAYAK digunakan';
  else
    status = 'TIDAK layak digunakan';
  end
  fprintf('   n = %d : rerata galat relatif %.6f %% -> %s\n', ...
          N(m), rerata_relatif(m), status);
end
fprintf('==================================================================\n');

% ---- 5) grafik: keempat fungsi dalam satu plot ----
warna = {'r', 'g', 'b'};
gaya  = {'--', '-.', ':'};

figure(1);
subplot(2,1,1);
plot(t, h_eksak, 'k-', 'LineWidth', 2.5); hold on;
for m = 1:numel(N)
  plot(t, h_aprox(m,:), [warna{m} gaya{m}], 'LineWidth', 1.5);
end
hold off; grid on;
xlabel('t'); ylabel('h(t)');
title('Perbandingan h(t): eksak vs pendekatan deret Taylor');
legend('f eksak', 'f n=23', 'f n=24', 'f n=25', 'Location', 'best');

% zoom supaya perbedaan antar metode terlihat jelas
subplot(2,1,2);
idx = t <= 0.5;
plot(t(idx), h_eksak(idx), 'k-', 'LineWidth', 2.5); hold on;
for m = 1:numel(N)
  plot(t(idx), h_aprox(m,idx), [warna{m} gaya{m}], 'LineWidth', 1.5);
end
hold off; grid on;
xlabel('t'); ylabel('h(t)');
title('Zoom pada 0 <= t <= 0.5');
legend('f eksak', 'f n=23', 'f n=24', 'f n=25', 'Location', 'best');
print(gcf, 'kasus1_grafik_perbandingan.png', '-dpng', '-S760,900');

% ---- 6) grafik galat mutlak terhadap t ----
figure(2);
for m = 1:numel(N)
  semilogy(t, abs(h_eksak - h_aprox(m,:)), [warna{m} gaya{m}], 'LineWidth', 1.5);
  hold on;
end
hold off; grid on;
xlabel('t'); ylabel('Galat mutlak');
title('Galat mutlak terhadap hasil eksak');
legend('n=23', 'n=24', 'n=25', 'Location', 'best');
print(gcf, 'kasus1_grafik_galat.png', '-dpng', '-r120');
