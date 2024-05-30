% Maximization Problems

function index = rouletteWheel( skor )
    r = rand * sum(skor);
    c = cumsum(skor);
    index = find( r <= c, 1, 'first');
end