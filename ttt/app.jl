using Gtk

include("main.jl")

win = GtkWindow("Tic-Tac-Toe", 400, 400)

newButt = GtkButton("New Game")

human = 1
player = 1

g = GtkGrid()
newButt = GtkButton("New Game")
g[1:3,1] = newButt

signal_connect(newButt, "clicked") do widget
    for c=1:3
        for r=2:4
            set_gtk_property!(g[c,r],:label," ")
        end
    end
    global player = 1
    global human = 3 - human
    if human==2
        _,best = findBest(position0)
        (r,c) = findfirst(map((x,y) -> x!=y,position0.pos,best.pos)).I
        set_gtk_property!(g[c,r+1],:label,"X")
        player = 2
    end
end

for c=1:3
    for r=2:4
        b = GtkButton(" ")
        g[c,r] = b
        signal_connect(b, "clicked") do widget
            oldsym = get_gtk_property(b,:label,String)
            if oldsym==" "
                global player
                global human
                sym = player==1 ? "X" : "O"
                set_gtk_property!(b,:label,sym)
                player = 3 - player
                mat,n = encodePos()
                p = position(mat,player,0,n)
                p = setStaticVal(p)
                if p.val!=0
                    postmessage("You win!")
                elseif p.n>9
                    postmessage("A Draw")
                else
                    _,best = findBest(p)
                    if best.val!=0
                        postmessage("I win!")
                    elseif p.n>9
                        postmessage("A Draw")
                    end
                    (r,c) = findfirst(map((x,y) -> x!=y,p.pos,best.pos)).I
                    sym = player==1 ? "X" : "O"
                    set_gtk_property!(g[c,r+1],:label,sym)
                    player = 3 - player
                end
            end
        end
    end
end
messageBox = GtkLabel("messages...")
g[1:3,5:7] = messageBox
push!(win, g)

showall(win)


function encodePos()::Tuple{Matrix{Int8},Int8}
    mat = zeros(Int8,3,3)
    n = 1
    for c = 1:3
        for r = 1:3
            sym = get_gtk_property(g[c,r+1],:label,String)
            if sym!=" "
                n += 1
                val = sym=="X" ? 1 : 2
                mat[r,c] = val
            end
        end
    end
    mat,n
end

function postmessage(mess::String)
    set_gtk_property!(messageBox,:label,mess)
end
