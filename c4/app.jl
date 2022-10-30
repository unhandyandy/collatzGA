using Gtk,JLD2,FileIO

include("engine.jl")

win = GtkWindow("Connect4", 400, 300)

newButt = GtkButton("New Game")

human = 1
node = makeStartNode()
depth = 10

g = GtkGrid()
newButt = GtkButton("New Game")
g[1:3,1] = newButt

undoButt = GtkButton("Undo")
g[5:7,1] = undoButt

@JLD2.load "memory.jld2"
@JLD2.load "coefficients.jld2"

function playerToPiece(p)
    p==1 ? "X" : "O"
end

signal_connect(newButt, "clicked") do widget
    ## post game stuff
    ## switch following lines to alternate X and O
    global human = 3 - human
    # global human
    ## in the following line human = comp in previous game
    postMortem(node, human)
    updateCoefficients()
    ## new game stuff
    postmessage(string("New game: you are ",playerToPiece(human)))
    global node = makeStartNode()
    for c=1:7
        for r=2:7
            set_gtk_property!(g[c,r],:label," ")
        end
    end
    if human==2
        compMove()
    end
end

function updateBoard()
    mat = node.position.board
    for c = 1:7
        for r = 1:6
            val = mat[r,c]
            sym = val==0 ? " " : val==1 ? "X" : "O"
            set_gtk_property!(g[c,r+1],:label,sym)
        end
    end
end

signal_connect(undoButt, "clicked") do widget
    global node
    postmessage("Undoing")
    display(node.parent)
    node = node.parent.parent
    updateBoard()
end

function replaceNode(newn::Node)
    global node
    newn.parent = node
    # println("parent set")
    node = newn
end

for c=1:7
    for r=2:7
        b = GtkButton(" ")
        g[c,r] = b
        signal_connect(b, "clicked") do widget
            # println(r,c)
            if isLegal(r-1,c)
                global human,node
                sym = node.position.mover==1 ? "X" : "O"
                set_gtk_property!(b,:label,sym)
                mat = encodePos()
                p = Position(mat,toggleMover(node.position.mover))
                newnode = Node(p,node.n+1,node.hist)
                replaceNode(newnode)
                display(node)
                # println("replaced")
                if !checkEnd(node)
                    compMove()
                end
            else
                postmessage("Illegal move")
            end
        end
    end
end
messageBox = GtkLabel("messages...")
g[1:3,9:10] = messageBox
push!(win, g)

showall(win)

function encodePos()
    mat = zeros(Int8,6,7)
    for c = 1:7
        for r = 1:6
            sym = get_gtk_property(g[c,r+1],:label,String)
            if sym!=" "
                val = sym=="X" ? 1 : 2
                mat[r,c] = val
            end
        end
    end
    mat
end

function postmessage(mess::String)
    set_gtk_property!(messageBox,:label,mess)
end

function checkEnd(n::Node)
    # display(n.position)
    if n.position.term
        if n.position.status==-1
            who = human==node.position.mover ? "I" : "You"
            postmessage(string(who," win!"))
        else
            postmessage("A Draw")
        end
        @JLD2.save "memory.jld2" memory
        @JLD2.save "coefficients.jld2" coefficients
        return true
    else
        return false
    end
end

function isLegal(r,c)
    return findHeight(node.position.board,c)==r && !node.position.term
end

function compMove()
    global node
    oldbd = node.position.board
    plyr = node.position.mover
    postmessage("Thinking...")
    node = chooseMove(node,depth)
    postmessage("...done.")
    display(node)
    (nr,nc) = findfirst(map((x,y) -> x!=y,oldbd,node.position.board)).I
    # println(nr,nc)
    sym = plyr==1 ? "X" : "O"
    set_gtk_property!(g[nc,nr+1],:label,sym)
    checkEnd(node)
end
