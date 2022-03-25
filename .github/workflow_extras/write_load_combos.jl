using Combinatorics
open("load_combos.sh", "w") do f
    for (a,b,c,d) in permutations(ARGS)
        write(f, "julia --project -e 'using $a; using $b; using $c; using $d' && \\\n")
    end
    write(f, "echo -e \\\\nDone!\n")
end